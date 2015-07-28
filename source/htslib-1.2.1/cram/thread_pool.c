/*
Copyright (c) 2013 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <signal.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <chrono>
#else
#include <sys/time.h>
#endif
#include <assert.h>

#include "cram/thread_pool.h"

//#define DEBUG
//#define DEBUG_TIME

#define IN_ORDER

/* Core of t_pool_next_result() */
static t_pool_result *t_pool_next_result_locked(t_results_queue *q) {
	t_pool_result *r, *last;

	for (last = NULL, r = q->result_head; r; last = r, r = r->next) {
		if (r->serial == q->next_serial)
			break;
	}

	if (r) {
		if (q->result_head == r)
			q->result_head = r->next;
		else
			last->next = r->next;

		if (q->result_tail == r)
			q->result_tail = last;

		if (!q->result_head)
			q->result_tail = NULL;

		q->next_serial++;
		q->queue_len--;
}

	return r;
}

/*
* Frees a result 'r' and if free_data is true also frees
* the internal r->data result too.
*/
void t_pool_delete_result(t_pool_result *r, int free_data) {
	if (!r)
		return;

	if (free_data && r->data)
		free(r->data);

	free(r);
}


#if !defined(_WIN32) && defined(USE_PTHREAD) 

#ifdef DEBUG
static int worker_id(t_pool *p) {
	int i;
	pthread_t s = pthread_self();
	for (i = 0; i < p->tsize; i++) {
		if (pthread_equal(s, p->t[i].tid))
			return i;
	}
	return -1;
}
#endif

/* ----------------------------------------------------------------------------
* A queue to hold results from the thread pool.
*
* Each thread pool may have jobs of multiple types being queued up and
* interleaved, so we allow several results queue per pool.
*
* The jobs themselves are expected to push their results onto their
* appropriate results queue.
*/

/*
* Adds a result to the end of the result queue.
*
* Returns 0 on success;
*        -1 on failure
*/
static int t_pool_add_result(t_pool_job *j, void *data) {
	t_results_queue *q = j->q;
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "%d: Adding resulting to queue %p, serial %d\n",
		worker_id(j->p), q, j->serial);
#endif

	/* No results queue is fine if we don't want any results back */
	if (!q)
		return 0;

	if (!(r = malloc(sizeof(*r))))
		return -1;

	r->next = NULL;
	r->data = data;
	r->serial = j->serial;

	pthread_mutex_lock(&q->result_m);
	if (q->result_tail) {
		q->result_tail->next = r;
		q->result_tail = r;
	}
	else {
		q->result_head = q->result_tail = r;
	}
	q->queue_len++;
	q->pending--;

#ifdef DEBUG
	fprintf(stderr, "%d: Broadcasting result_avail (id %d)\n",
		worker_id(j->p), r->serial);
#endif
	pthread_cond_signal(&q->result_avail_c);
#ifdef DEBUG
	fprintf(stderr, "%d: Broadcast complete\n", worker_id(j->p));
#endif

	pthread_mutex_unlock(&q->result_m);

	return 0;
}

/*
* Pulls a result off the head of the result queue. Caller should
* free it (and any internals as appropriate) after use. This doesn't
* wait for a result to be present.
*
* Results will be returned in strict order.
*
* Returns t_pool_result pointer if a result is ready.
*         NULL if not.
*/
t_pool_result *t_pool_next_result(t_results_queue *q) {
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "Requesting next result on queue %p\n", q);
#endif

	pthread_mutex_lock(&q->result_m);
	r = t_pool_next_result_locked(q);
	pthread_mutex_unlock(&q->result_m);

#ifdef DEBUG
	fprintf(stderr, "(q=%p) Found %p\n", q, r);
#endif

	return r;
}

t_pool_result *t_pool_next_result_wait(t_results_queue *q) {
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "Waiting for result %d...\n", q->next_serial);
#endif

	pthread_mutex_lock(&q->result_m);
	while (!(r = t_pool_next_result_locked(q))) {
		/* Possible race here now avoided via _locked() call, but incase... */
		struct timeval now;
		struct timespec timeout;

		gettimeofday(&now, NULL);
		timeout.tv_sec = now.tv_sec + 10;
		timeout.tv_nsec = now.tv_usec * 1000;

		pthread_cond_timedwait(&q->result_avail_c, &q->result_m, &timeout);
	}
	pthread_mutex_unlock(&q->result_m);

	return r;
}

/*
* Returns true if there are no items on the finished results queue and
* also none still pending.
*/
int t_pool_results_queue_empty(t_results_queue *q) {
	int empty;

	pthread_mutex_lock(&q->result_m);
	empty = q->queue_len == 0 && q->pending == 0;
	pthread_mutex_unlock(&q->result_m);

	return empty;
}


/*
* Returns the number of completed jobs on the results queue.
*/
int t_pool_results_queue_len(t_results_queue *q) {
	int len;

	pthread_mutex_lock(&q->result_m);
	len = q->queue_len;
	pthread_mutex_unlock(&q->result_m);

	return len;
}

int t_pool_results_queue_sz(t_results_queue *q) {
	int len;

	pthread_mutex_lock(&q->result_m);
	len = q->queue_len + q->pending;
	pthread_mutex_unlock(&q->result_m);

	return len;
}

/*
* Initialises a results queue.
*
* Results queue pointer on success;
*         NULL on failure
*/
t_results_queue *t_results_queue_init(void) {
	t_results_queue *q = malloc(sizeof(*q));

	pthread_mutex_init(&q->result_m, NULL);
	pthread_cond_init(&q->result_avail_c, NULL);

	q->result_head = NULL;
	q->result_tail = NULL;
	q->next_serial = 0;
	q->curr_serial = 0;
	q->queue_len = 0;
	q->pending = 0;

	return q;
}

/* Deallocates memory for a results queue */
void t_results_queue_destroy(t_results_queue *q) {
#ifdef DEBUG
	fprintf(stderr, "Destroying results queue %p\n", q);
#endif

	if (!q)
		return;

	pthread_mutex_destroy(&q->result_m);
	pthread_cond_destroy(&q->result_avail_c);

	memset(q, 0xbb, sizeof(*q));
	free(q);

#ifdef DEBUG
	fprintf(stderr, "Destroyed results queue %p\n", q);
#endif
}

/* ----------------------------------------------------------------------------
* The thread pool.
*/

#define TDIFF(t2,t1) ((t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec)

/*
* A worker thread.
*
* Each thread waits for the pool to be non-empty.
* As soon as this applies, one of them succeeds in getting the lock
* and then executes the job.
*/
static void *t_pool_worker(void *arg) {
	t_pool_worker_t *w = (t_pool_worker_t *)arg;
	t_pool *p = w->p;
	t_pool_job *j;
#ifdef DEBUG_TIME
	struct timeval t1, t2, t3;
#endif

	for (;;) {
		// Pop an item off the pool queue
#ifdef DEBUG_TIME
		gettimeofday(&t1, NULL);
#endif

		pthread_mutex_lock(&p->pool_m);

#ifdef DEBUG_TIME
		gettimeofday(&t2, NULL);
		p->wait_time += TDIFF(t2, t1);
		w->wait_time += TDIFF(t2, t1);
#endif

		// If there is something on the job list and a higher priority
		// thread waiting, let it handle this instead.
		//	while (p->head && p->t_stack_top != -1 && p->t_stack_top < w->idx) {
		//	    pthread_mutex_unlock(&p->pool_m);
		//	    pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
		//	    pthread_mutex_lock(&p->pool_m);
		//	}

		while (!p->head && !p->shutdown) {
			p->nwaiting++;

			if (p->njobs == 0)
				pthread_cond_signal(&p->empty_c);
#ifdef DEBUG_TIME
			gettimeofday(&t2, NULL);
#endif

#ifdef IN_ORDER
			// Push this thread to the top of the waiting stack
			if (p->t_stack_top == -1 || p->t_stack_top > w->idx)
				p->t_stack_top = w->idx;

			p->t_stack[w->idx] = 1;
			pthread_cond_wait(&w->pending_c, &p->pool_m);
			p->t_stack[w->idx] = 0;

			/* Find new t_stack_top */
			{
				int i;
				p->t_stack_top = -1;
				for (i = 0; i < p->tsize; i++) {
					if (p->t_stack[i]) {
						p->t_stack_top = i;
						break;
					}
				}
			}
#else
			pthread_cond_wait(&p->pending_c, &p->pool_m);
#endif

#ifdef DEBUG_TIME
			gettimeofday(&t3, NULL);
			p->wait_time += TDIFF(t3, t2);
			w->wait_time += TDIFF(t3, t2);
#endif
			p->nwaiting--;
		}

		if (p->shutdown) {
#ifdef DEBUG_TIME
			p->total_time += TDIFF(t3, t1);
#endif
#ifdef DEBUG
			fprintf(stderr, "%d: Shutting down\n", worker_id(p));
#endif
			pthread_mutex_unlock(&p->pool_m);
			pthread_exit(NULL);
		}

		j = p->head;
		if (!(p->head = j->next))
			p->tail = NULL;

		if (p->njobs-- >= p->qsize)
			pthread_cond_signal(&p->full_c);

		if (p->njobs == 0)
			pthread_cond_signal(&p->empty_c);

		pthread_mutex_unlock(&p->pool_m);

		// We have job 'j' - now execute it.
		t_pool_add_result(j, j->func(j->arg));
#ifdef DEBUG_TIME
		pthread_mutex_lock(&p->pool_m);
		gettimeofday(&t3, NULL);
		p->total_time += TDIFF(t3, t1);
		pthread_mutex_unlock(&p->pool_m);
#endif
		memset(j, 0xbb, sizeof(*j));
		free(j);
	}

	return NULL;
}

/*
* Creates a worker pool of length qsize with tsize worker threads.
*
* Returns pool pointer on success;
*         NULL on failure
*/
t_pool *t_pool_init(int qsize, int tsize) {
	int i;
	t_pool *p = malloc(sizeof(*p));
	p->qsize = qsize;
	p->tsize = tsize;
	p->njobs = 0;
	p->nwaiting = 0;
	p->shutdown = 0;
	p->head = p->tail = NULL;
	p->t_stack = NULL;
#ifdef DEBUG_TIME
	p->total_time = p->wait_time = 0;
#endif

	p->t = malloc(tsize * sizeof(p->t[0]));

	pthread_mutex_init(&p->pool_m, NULL);
	pthread_cond_init(&p->empty_c, NULL);
	pthread_cond_init(&p->full_c, NULL);

	pthread_mutex_lock(&p->pool_m);

#ifdef IN_ORDER
	if (!(p->t_stack = malloc(tsize * sizeof(*p->t_stack))))
		return NULL;
	p->t_stack_top = -1;

	for (i = 0; i < tsize; i++) {
		t_pool_worker_t *w = &p->t[i];
		p->t_stack[i] = 0;
		w->p = p;
		w->idx = i;
		w->wait_time = 0;
		pthread_cond_init(&w->pending_c, NULL);
		if (0 != pthread_create(&w->tid, NULL, t_pool_worker, w))
			return NULL;
	}
#else
	pthread_cond_init(&p->pending_c, NULL);

	for (i = 0; i < tsize; i++) {
		t_pool_worker_t *w = &p->t[i];
		w->p = p;
		w->idx = i;
		pthread_cond_init(&w->pending_c, NULL);
		if (0 != pthread_create(&w->tid, NULL, t_pool_worker, w))
			return NULL;
	}
#endif

	pthread_mutex_unlock(&p->pool_m);

	return p;
}

/*
* Adds an item to the work pool.
*
* FIXME: Maybe return 1,0,-1 and distinguish between job dispathed vs
* result returned. Ie rather than blocking on full queue we're permitted
* to return early on "result available" event too.
* Caller would then have a while loop around t_pool_dispatch.
* Or, return -1 and set errno to EAGAIN to indicate job not yet submitted.
*
* Returns 0 on success
*        -1 on failure
*/
int t_pool_dispatch(t_pool *p, t_results_queue *q,
	void *(*func)(void *arg), void *arg) {
	t_pool_job *j = malloc(sizeof(*j));

	if (!j)
		return -1;
	j->func = func;
	j->arg = arg;
	j->next = NULL;
	j->p = p;
	j->q = q;
	if (q) {
		pthread_mutex_lock(&q->result_m);
		j->serial = q->curr_serial++;
		q->pending++;
		pthread_mutex_unlock(&q->result_m);
	}
	else {
		j->serial = 0;
	}

#ifdef DEBUG
	fprintf(stderr, "Dispatching job %p for queue %p, serial %d\n", j, q, j->serial);
#endif

	pthread_mutex_lock(&p->pool_m);

	// Check if queue is full
	while (p->njobs >= p->qsize)
		pthread_cond_wait(&p->full_c, &p->pool_m);

	p->njobs++;

	if (p->tail) {
		p->tail->next = j;
		p->tail = j;
	}
	else {
		p->head = p->tail = j;
	}

	// Let a worker know we have data.
#ifdef IN_ORDER
	if (p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting)
		pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
#else
	pthread_cond_signal(&p->pending_c);
#endif
	pthread_mutex_unlock(&p->pool_m);

#ifdef DEBUG
	fprintf(stderr, "Dispatched (serial %d)\n", j->serial);
#endif

	return 0;
}

/*
* As above but optional non-block flag.
*
* nonblock  0 => block if input queue is full
* nonblock +1 => don't block if input queue is full, but do not add task
* nonblock -1 => add task regardless of whether queue is full (over-size)
*/
int t_pool_dispatch2(t_pool *p, t_results_queue *q,
	void *(*func)(void *arg), void *arg, int nonblock) {
	t_pool_job *j;

#ifdef DEBUG
	fprintf(stderr, "Dispatching job for queue %p, serial %d\n", q, q->curr_serial);
#endif

	pthread_mutex_lock(&p->pool_m);

	if (p->njobs >= p->qsize && nonblock == 1) {
		pthread_mutex_unlock(&p->pool_m);
		errno = EAGAIN;
		return -1;
	}

	if (!(j = malloc(sizeof(*j))))
		return -1;
	j->func = func;
	j->arg = arg;
	j->next = NULL;
	j->p = p;
	j->q = q;
	if (q) {
		pthread_mutex_lock(&q->result_m);
		j->serial = q->curr_serial;
		pthread_mutex_unlock(&q->result_m);
	}
	else {
		j->serial = 0;
	}

	if (q) {
		pthread_mutex_lock(&q->result_m);
		q->curr_serial++;
		q->pending++;
		pthread_mutex_unlock(&q->result_m);
	}

	// Check if queue is full
	if (nonblock == 0)
		while (p->njobs >= p->qsize)
			pthread_cond_wait(&p->full_c, &p->pool_m);

	p->njobs++;

	//    if (q->curr_serial % 100 == 0)
	//	fprintf(stderr, "p->njobs = %d    p->qsize = %d\n", p->njobs, p->qsize);

	if (p->tail) {
		p->tail->next = j;
		p->tail = j;
	}
	else {
		p->head = p->tail = j;
	}

#ifdef DEBUG
	fprintf(stderr, "Dispatched (serial %d)\n", j->serial);
#endif

	// Let a worker know we have data.
#ifdef IN_ORDER
	// Keep incoming queue at 1 per running thread, so there is always
	// something waiting when they end their current task.  If we go above
	// this signal to start more threads (if available). This has the effect
	// of concentrating jobs to fewer cores when we are I/O bound, which in
	// turn benefits systems with auto CPU frequency scaling.
	if (p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting)
		pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
#else
	pthread_cond_signal(&p->pending_c);
#endif

	pthread_mutex_unlock(&p->pool_m);

	return 0;
}

/*
* Flushes the pool, but doesn't exit. This simply drains the queue and
* ensures all worker threads have finished their current task.
*
* Returns 0 on success;
*        -1 on failure
*/
int t_pool_flush(t_pool *p) {
	int i;

#ifdef DEBUG
	fprintf(stderr, "Flushing pool %p\n", p);
#endif

	// Drains the queue
	pthread_mutex_lock(&p->pool_m);

	// Wake up everything for the final sprint!
	for (i = 0; i < p->tsize; i++)
		if (p->t_stack[i])
			pthread_cond_signal(&p->t[i].pending_c);

	while (p->njobs || p->nwaiting != p->tsize)
		pthread_cond_wait(&p->empty_c, &p->pool_m);

	pthread_mutex_unlock(&p->pool_m);

#ifdef DEBUG
	fprintf(stderr, "Flushed complete for pool %p, njobs=%d, nwaiting=%d\n",
		p, p->njobs, p->nwaiting);
#endif

	return 0;
}

/*
* Destroys a thread pool. If 'kill' is true the threads are terminated now,
* otherwise they are joined into the main thread so they will finish their
* current work load.
*
* Use t_pool_destroy(p,0) after a t_pool_flush(p) on a normal shutdown or
* t_pool_destroy(p,1) to quickly exit after a fatal error.
*/
void t_pool_destroy(t_pool *p, int kill) {
	int i;

#ifdef DEBUG
	fprintf(stderr, "Destroying pool %p, kill=%d\n", p, kill);
#endif

	/* Send shutdown message to worker threads */
	if (!kill) {
		pthread_mutex_lock(&p->pool_m);
		p->shutdown = 1;

#ifdef DEBUG
		fprintf(stderr, "Sending shutdown request\n");
#endif

#ifdef IN_ORDER
		for (i = 0; i < p->tsize; i++)
			pthread_cond_signal(&p->t[i].pending_c);
#else
		pthread_cond_broadcast(&p->pending_c);
#endif
		pthread_mutex_unlock(&p->pool_m);

#ifdef DEBUG
		fprintf(stderr, "Shutdown complete\n");
#endif
		for (i = 0; i < p->tsize; i++)
			pthread_join(p->t[i].tid, NULL);
	}
	else {
		for (i = 0; i < p->tsize; i++)
			pthread_kill(p->t[i].tid, SIGINT);
	}

	pthread_mutex_destroy(&p->pool_m);
	pthread_cond_destroy(&p->empty_c);
	pthread_cond_destroy(&p->full_c);
#ifdef IN_ORDER
	for (i = 0; i < p->tsize; i++)
		pthread_cond_destroy(&p->t[i].pending_c);
#else
	pthread_cond_destroy(&p->pending_c);
#endif

#ifdef DEBUG_TIME
	fprintf(stderr, "Total time=%f\n", p->total_time / 1000000.0);
	fprintf(stderr, "Wait  time=%f\n", p->wait_time / 1000000.0);
	fprintf(stderr, "%d%% utilisation\n",
		(int)(100 - ((100.0 * p->wait_time) / p->total_time + 0.5)));
	for (i = 0; i < p->tsize; i++)
		fprintf(stderr, "%d: Wait time=%f\n", i,
		p->t[i].wait_time / 1000000.0);
#endif

	if (p->t_stack)
		free(p->t_stack);

	free(p->t);
	free(p);

#ifdef DEBUG
	fprintf(stderr, "Destroyed pool %p\n", p);
#endif
}

#else // ~ !defined(_WIN32) && defined(USE_PTHREAD) 

#ifdef DEBUG
static int worker_id(t_pool *p) {
	int i;
	pthread_t s = pthread_self();
	for (i = 0; i < p->tsize; i++) {
		if (pthread_equal(s, p->t[i].tid))
			return i;
	}
	return -1;
}
#endif

/* ----------------------------------------------------------------------------
* A queue to hold results from the thread pool.
*
* Each thread pool may have jobs of multiple types being queued up and
* interleaved, so we allow several results queue per pool.
*
* The jobs themselves are expected to push their results onto their
* appropriate results queue.
*/

/*
* Adds a result to the end of the result queue.
*
* Returns 0 on success;
*        -1 on failure
*/
static int t_pool_add_result(t_pool_job *j, void *data) {
	t_results_queue *q = j->q;
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "%d: Adding resulting to queue %p, serial %d\n",
		worker_id(j->p), q, j->serial);
#endif

	/* No results queue is fine if we don't want any results back */
	if (!q)
		return 0;

	if (!(r = (t_pool_result*)malloc(sizeof(*r))))
		return -1;

	r->next = NULL;
	r->data = data;
	r->serial = j->serial;
	q->result_m.lock();
	if (q->result_tail) {
		q->result_tail->next = r;
		q->result_tail = r;
	}
	else {
		q->result_head = q->result_tail = r;
	}
	q->queue_len++;
	q->pending--;

#ifdef DEBUG
	fprintf(stderr, "%d: Broadcasting result_avail (id %d)\n",
		worker_id(j->p), r->serial);
#endif
	q->result_avail_c.notify_one();
#ifdef DEBUG
	fprintf(stderr, "%d: Broadcast complete\n", worker_id(j->p));
#endif
	q->result_m.unlock();
	return 0;
}

/*
* Pulls a result off the head of the result queue. Caller should
* free it (and any internals as appropriate) after use. This doesn't
* wait for a result to be present.
*
* Results will be returned in strict order.
*
* Returns t_pool_result pointer if a result is ready.
*         NULL if not.
*/
t_pool_result *t_pool_next_result(t_results_queue *q) {
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "Requesting next result on queue %p\n", q);
#endif

	q->result_m.lock();
	r = t_pool_next_result_locked(q);
	q->result_m.unlock();

#ifdef DEBUG
	fprintf(stderr, "(q=%p) Found %p\n", q, r);
#endif

	return r;
}


t_pool_result *t_pool_next_result_wait(t_results_queue *q) {
	t_pool_result *r;

#ifdef DEBUG
	fprintf(stderr, "Waiting for result %d...\n", q->next_serial);
#endif

	q->result_m.lock();
	while (!(r = t_pool_next_result_locked(q))) {
		/* Possible race here now avoided via _locked() call, but incase... */

		std::unique_lock<std::mutex> lk(q->result_m);
		auto now = std::chrono::system_clock::now();
		std::chrono::seconds sec(10); 
		
		q->result_avail_c.wait_until(lk, now + sec); 

	}
	q->result_m.unlock();
	return r;
}

/*
* Returns true if there are no items on the finished results queue and
* also none still pending.
*/
int t_pool_results_queue_empty(t_results_queue *q) {
	int empty;

	q->result_m.lock();
	empty = q->queue_len == 0 && q->pending == 0;
	q->result_m.unlock();

	return empty;
}

/*
* Returns the number of completed jobs on the results queue.
*/
int t_pool_results_queue_len(t_results_queue *q) {
	int len;

	q->result_m.lock();
	len = q->queue_len;
	q->result_m.unlock();

	return len;
}

int t_pool_results_queue_sz(t_results_queue *q) {
	int len;

	q->result_m.lock();
	len = q->queue_len + q->pending;
	q->result_m.unlock();

	return len;
}

#define TDIFF(t2,t1) ((t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec)

/*
* Initialises a results queue.
*
* Results queue pointer on success;
*         NULL on failure
*/
t_results_queue *t_results_queue_init(void) {
	t_results_queue *q = (t_results_queue*)malloc(sizeof(*q));

	q->result_head = NULL;
	q->result_tail = NULL;
	q->next_serial = 0;
	q->curr_serial = 0;
	q->queue_len = 0;
	q->pending = 0;

	return q;
}

/* Deallocates memory for a results queue */
void t_results_queue_destroy(t_results_queue *q) {
#ifdef DEBUG
	fprintf(stderr, "Destroying results queue %p\n", q);
#endif

	if (!q)
		return;

	memset(q, 0xbb, sizeof(*q));
	free(q);

#ifdef DEBUG
	fprintf(stderr, "Destroyed results queue %p\n", q);
#endif
}

/* ----------------------------------------------------------------------------
* The thread pool.
*/

#define TDIFF(t2,t1) ((t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec)

/*
* A worker thread.
*
* Each thread waits for the pool to be non-empty.
* As soon as this applies, one of them succeeds in getting the lock
* and then executes the job.
*/
static void *t_pool_worker(void *arg) {
	t_pool_worker_t *w = (t_pool_worker_t *)arg;
	t_pool *p = w->p;
	t_pool_job *j;
#ifdef DEBUG_TIME
	struct timeval t1, t2, t3;
#endif

	for (;;) {
		// Pop an item off the pool queue
#ifdef DEBUG_TIME
		gettimeofday(&t1, NULL);
#endif

		p->pool_m.lock();

#ifdef DEBUG_TIME
		gettimeofday(&t2, NULL);
		p->wait_time += TDIFF(t2, t1);
		w->wait_time += TDIFF(t2, t1);
#endif

		// If there is something on the job list and a higher priority
		// thread waiting, let it handle this instead.
		//	while (p->head && p->t_stack_top != -1 && p->t_stack_top < w->idx) {
		//	    pthread_mutex_unlock(&p->pool_m);
		//	    pthread_cond_signal(&p->t[p->t_stack_top].pending_c);
		//	    pthread_mutex_lock(&p->pool_m);
		//	}

		while (!p->head && !p->shutdown) {
			p->nwaiting++;

			if (p->njobs == 0)
				p->empty_c.notify_one();
#ifdef DEBUG_TIME
			gettimeofday(&t2, NULL);
#endif

#ifdef IN_ORDER
			// Push this thread to the top of the waiting stack
			if (p->t_stack_top == -1 || p->t_stack_top > w->idx)
				p->t_stack_top = w->idx;

			p->t_stack[w->idx] = 1;
			std::unique_lock<std::mutex> lk(p->pool_m); 
			w->pending_c.wait(lk); 
			p->t_stack[w->idx] = 0;

			/* Find new t_stack_top */
			{
				int i;
				p->t_stack_top = -1;
				for (i = 0; i < p->tsize; i++) {
					if (p->t_stack[i]) {
						p->t_stack_top = i;
						break;
					}
				}
			}
#else
			std::unique_lock<std::mutex> lk(p->pool_m); 
			p->pending_c.wait(lk);
#endif

#ifdef DEBUG_TIME
			gettimeofday(&t3, NULL);
			p->wait_time += TDIFF(t3, t2);
			w->wait_time += TDIFF(t3, t2);
#endif
			p->nwaiting--;
		}

		if (p->shutdown) {
#ifdef DEBUG_TIME
			p->total_time += TDIFF(t3, t1);
#endif
#ifdef DEBUG
			fprintf(stderr, "%d: Shutting down\n", worker_id(p));
#endif
			p->pool_m.unlock();
			// return function to exit thread, no specific method to exit std::thread
			return (NULL);
		}

		j = p->head;
		if (!(p->head = j->next))
			p->tail = NULL;

		if (p->njobs-- >= p->qsize)
			p->full_c.notify_one();

		if (p->njobs == 0)
			p->empty_c.notify_one();

		p->pool_m.unlock();

		// We have job 'j' - now execute it.
		t_pool_add_result(j, j->func(j->arg));
#ifdef DEBUG_TIME
		p->pool_m.lock();
		gettimeofday(&t3, NULL);
		p->total_time += TDIFF(t3, t1);
		p->pool_m.unlock();
#endif
		memset(j, 0xbb, sizeof(*j));
		free(j);
	}

	return NULL;
}

/*
* Creates a worker pool of length qsize with tsize worker threads.
*
* Returns pool pointer on success;
*         NULL on failure
*/
t_pool *t_pool_init(int qsize, int tsize) {
	int i;
	t_pool *p = (t_pool*)malloc(sizeof(*p));
	p->qsize = qsize;
	p->tsize = tsize;
	p->njobs = 0;
	p->nwaiting = 0;
	p->shutdown = 0;
	p->head = p->tail = NULL;
	p->t_stack = NULL;
#ifdef DEBUG_TIME
	p->total_time = p->wait_time = 0;
#endif

	p->t = (t_pool_worker_t*)malloc(tsize * sizeof(p->t[0]));

	p->pool_m.lock();

#ifdef IN_ORDER
	if (!(p->t_stack = (int*)malloc(tsize * sizeof(*p->t_stack))))
		return NULL;
	p->t_stack_top = -1;

	for (i = 0; i < tsize; i++) {
		t_pool_worker_t *w = &p->t[i];
		p->t_stack[i] = 0;
		w->p = p;
		w->idx = i;
		w->wait_time = 0;

		p->threads.push_back(std::thread(t_pool_worker, w)); 
	}
#else

	for (i = 0; i < tsize; i++) {
		t_pool_worker_t *w = &p->t[i];
		w->p = p;
		w->idx = i;

		p->threads.push_back(std::thread(t_pool_worker, w)); 
	}
#endif

	p->pool_m.unlock();

	return p;
}

/*
* Adds an item to the work pool.
*
* FIXME: Maybe return 1,0,-1 and distinguish between job dispathed vs
* result returned. Ie rather than blocking on full queue we're permitted
* to return early on "result available" event too.
* Caller would then have a while loop around t_pool_dispatch.
* Or, return -1 and set errno to EAGAIN to indicate job not yet submitted.
*
* Returns 0 on success
*        -1 on failure
*/
int t_pool_dispatch(t_pool *p, t_results_queue *q,
	void *(*func)(void *arg), void *arg) {
	t_pool_job *j = (t_pool_job*)malloc(sizeof(*j));

	if (!j)
		return -1;
	j->func = func;
	j->arg = arg;
	j->next = NULL;
	j->p = p;
	j->q = q;
	if (q) {
		q->result_m.lock();
		j->serial = q->curr_serial++;
		q->pending++;
		q->result_m.unlock();
	}
	else {
		j->serial = 0;
	}

#ifdef DEBUG
	fprintf(stderr, "Dispatching job %p for queue %p, serial %d\n", j, q, j->serial);
#endif

	p->pool_m.lock();

	// Check if queue is full
	while (p->njobs >= p->qsize)
	{
		std::unique_lock<std::mutex> lk(p->pool_m);
		p->full_c.wait(lk);
	}
	p->njobs++;

	if (p->tail) {
		p->tail->next = j;
		p->tail = j;
	}
	else {
		p->head = p->tail = j;
	}

	// Let a worker know we have data.
#ifdef IN_ORDER
	if (p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting)
		p->t[p->t_stack_top].pending_c.notify_one();
#else
	p->pending_c.notify_one();
#endif
	p->pool_m.unlock();

#ifdef DEBUG
	fprintf(stderr, "Dispatched (serial %d)\n", j->serial);
#endif

	return 0;
}

/*
* As above but optional non-block flag.
*
* nonblock  0 => block if input queue is full
* nonblock +1 => don't block if input queue is full, but do not add task
* nonblock -1 => add task regardless of whether queue is full (over-size)
*/
int t_pool_dispatch2(t_pool *p, t_results_queue *q,
	void *(*func)(void *arg), void *arg, int nonblock) {
	t_pool_job *j;

#ifdef DEBUG
	fprintf(stderr, "Dispatching job for queue %p, serial %d\n", q, q->curr_serial);
#endif

	p->pool_m.lock();

	if (p->njobs >= p->qsize && nonblock == 1) {
		p->pool_m.unlock();
		errno = EAGAIN;
		return -1;
	}

	if (!(j = (t_pool_job*)malloc(sizeof(*j))))
		return -1;
	j->func = func;
	j->arg = arg;
	j->next = NULL;
	j->p = p;
	j->q = q;
	if (q) {
		q->result_m.lock();
		j->serial = q->curr_serial;
		q->result_m.unlock();
	}
	else {
		j->serial = 0;
	}

	if (q) {
		q->result_m.lock();
		q->curr_serial++;
		q->pending++;
		q->result_m.unlock();
	}

	// Check if queue is full
	if (nonblock == 0)
		while (p->njobs >= p->qsize)
		{
			std::unique_lock<std::mutex> lk(p->pool_m);
			p->full_c.wait(lk);
		}
	p->njobs++;

	//    if (q->curr_serial % 100 == 0)
	//	fprintf(stderr, "p->njobs = %d    p->qsize = %d\n", p->njobs, p->qsize);

	if (p->tail) {
		p->tail->next = j;
		p->tail = j;
	}
	else {
		p->head = p->tail = j;
	}

#ifdef DEBUG
	fprintf(stderr, "Dispatched (serial %d)\n", j->serial);
#endif

	// Let a worker know we have data.
#ifdef IN_ORDER
	// Keep incoming queue at 1 per running thread, so there is always
	// something waiting when they end their current task.  If we go above
	// this signal to start more threads (if available). This has the effect
	// of concentrating jobs to fewer cores when we are I/O bound, which in
	// turn benefits systems with auto CPU frequency scaling.
	if (p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting)
		p->t[p->t_stack_top].pending_c.notify_one();
#else
	p->pending_c.notify_one();
#endif

	p->pool_m.unlock();

	return 0;
}

/*
* Flushes the pool, but doesn't exit. This simply drains the queue and
* ensures all worker threads have finished their current task.
*
* Returns 0 on success;
*        -1 on failure
*/
int t_pool_flush(t_pool *p) {
	int i;

#ifdef DEBUG
	fprintf(stderr, "Flushing pool %p\n", p);
#endif

	// Drains the queue
	p->pool_m.lock();

	// Wake up everything for the final sprint!
	for (i = 0; i < p->tsize; i++)
		if (p->t_stack[i])
			p->t[i].pending_c.notify_one();

	while (p->njobs || p->nwaiting != p->tsize)
	{
		std::unique_lock<std::mutex> lk(p->pool_m);
		p->empty_c.wait(lk);
	}

	p->pool_m.unlock();

#ifdef DEBUG
	fprintf(stderr, "Flushed complete for pool %p, njobs=%d, nwaiting=%d\n",
		p, p->njobs, p->nwaiting);
#endif

	return 0;
}

/*
* Destroys a thread pool. If 'kill' is true the threads are terminated now,
* otherwise they are joined into the main thread so they will finish their
* current work load.
*
* Use t_pool_destroy(p,0) after a t_pool_flush(p) on a normal shutdown or
* t_pool_destroy(p,1) to quickly exit after a fatal error.
*/
void t_pool_destroy(t_pool *p, int kill) {
	int i;

#ifdef DEBUG
	fprintf(stderr, "Destroying pool %p, kill=%d\n", p, kill);
#endif

	/* Send shutdown message to worker threads */
	if (!kill) {
		p->pool_m.lock();
		p->shutdown = 1;

#ifdef DEBUG
		fprintf(stderr, "Sending shutdown request\n");
#endif

#ifdef IN_ORDER
		for (i = 0; i < p->tsize; i++)
			p->t[i].pending_c.notify_one();
#else
		p->pending_c.notify_all();
#endif
		p->pool_m.unlock();

#ifdef DEBUG
		fprintf(stderr, "Shutdown complete\n");
#endif
		for (auto& t : p->threads)
			t.join(); 
	}
	else {
		// TODO : how to kill std::thread ? 
		/*for (i = 0; i < p->tsize; i++)
			pthread_kill(p->t[i].tid, SIGINT);*/
	}

#ifdef DEBUG_TIME
	fprintf(stderr, "Total time=%f\n", p->total_time / 1000000.0);
	fprintf(stderr, "Wait  time=%f\n", p->wait_time / 1000000.0);
	fprintf(stderr, "%d%% utilisation\n",
		(int)(100 - ((100.0 * p->wait_time) / p->total_time + 0.5)));
	for (i = 0; i < p->tsize; i++)
		fprintf(stderr, "%d: Wait time=%f\n", i,
		p->t[i].wait_time / 1000000.0);
#endif

	if (p->t_stack)
		free(p->t_stack);

	free(p->t);
	free(p);

#ifdef DEBUG
	fprintf(stderr, "Destroyed pool %p\n", p);
#endif
}

#endif // ~ !defined(_WIN32) && defined(USE_PTHREAD) 


/*-----------------------------------------------------------------------------
* Test app.
*/

#ifdef TEST_MAIN

#include <stdio.h>
#include <math.h>

void *doit(void *arg) {
	int i, k, x = 0;
	int job = *(int *)arg;
	int *res;

	printf("Worker: execute job %d\n", job);

	usleep(random() % 1000000); // to coerce job completion out of order
	if (0) {
		for (k = 0; k < 100; k++) {
			for (i = 0; i < 100000; i++) {
				x++;
				x += x * sin(i);
				x += x * cos(x);
			}
		}
		x *= 100;
		x += job;
	}
	else {
		x = job*job;
	}

	printf("Worker: job %d terminating, x=%d\n", job, x);

	free(arg);

	res = malloc(sizeof(*res));
	*res = x;

	return res;
}

#define NTHREADS 8

int main(int argc, char **argv) {
	t_pool *p = t_pool_init(NTHREADS * 2, NTHREADS);
	t_results_queue *q = t_results_queue_init();
	int i;
	t_pool_result *r;

	// Dispatch jobs
	for (i = 0; i < 20; i++) {
		int *ip = malloc(sizeof(*ip));
		*ip = i;
		printf("Submitting %d\n", i);
		t_pool_dispatch(p, q, doit, ip);

		// Check for results
		if ((r = t_pool_next_result(q))) {
			printf("RESULT: %d\n", *(int *)r->data);
			t_pool_delete_result(r, 1);
		}
	}

	t_pool_flush(p);

	while ((r = t_pool_next_result(q))) {
		printf("RESULT: %d\n", *(int *)r->data);
		t_pool_delete_result(r, 1);
	}

	t_pool_destroy(p, 0);
	t_results_queue_destroy(q);

	return 0;
}
#endif