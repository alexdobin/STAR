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
#include <sys/time.h>

#include "cram/thread_pool.h"

//#define DEBUG
#define DEBUG_TIME

#ifdef DEBUG
static int worker_id(t_pool *p) {
    int i;
    pthread_t s = pthread_self();
    for (i = 0; i < p->tsize; i++) {
	if (pthread_equal(s, p->t[i]))
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
    } else {
	q->result_head = q->result_tail = r;
    }
    q->queue_len++;
    q->pending--;

#ifdef DEBUG
    fprintf(stderr, "%d: Broadcasting result_avail (id %d)\n",
	    worker_id(j->p), r->serial);
#endif
    pthread_cond_broadcast(&q->result_avail_c);
#ifdef DEBUG
    fprintf(stderr, "%d: Broadcast complete\n", worker_id(j->p));
#endif

    pthread_mutex_unlock(&q->result_m);

    return 0;
}

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
    q->queue_len   = 0;
    q->pending     = 0;

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
    t_pool *p = (t_pool *)arg;
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
	p->wait_time += TDIFF(t2,t1);
#endif

	p->nwaiting++;
	while (!p->head && !p->shutdown) {
	    if (p->njobs == 0)
		pthread_cond_signal(&p->empty_c);
#ifdef DEBUG_TIME
	    gettimeofday(&t2, NULL);
#endif

	    pthread_cond_wait(&p->pending_c, &p->pool_m);

#ifdef DEBUG_TIME
	    gettimeofday(&t3, NULL);
	    p->wait_time += TDIFF(t3,t2);
#endif
	}

	p->nwaiting--;

	if (p->shutdown) {
	    p->total_time += TDIFF(t3,t1);
#ifdef DEBUG
	    fprintf(stderr, "%d: Shutting down\n", worker_id(p));
#endif
	    pthread_mutex_unlock(&p->pool_m);
	    pthread_exit(NULL);
	}

	j = p->head;
	if (!(p->head = j->next))
	    p->tail = NULL;

	if (p->njobs-- == p->qsize)
	    pthread_cond_signal(&p->full_c);

	if (p->njobs == 0)
	    pthread_cond_signal(&p->empty_c);

	pthread_mutex_unlock(&p->pool_m);
	    
	// We have job 'j' - now execute it.
	t_pool_add_result(j, j->func(j->arg));	
#ifdef DEBUG_TIME
	pthread_mutex_lock(&p->pool_m);
	gettimeofday(&t3, NULL);
	p->total_time += TDIFF(t3,t1);
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
#ifdef DEBUG_TIME
    p->total_time = p->wait_time = 0;
#endif

    p->t = malloc(tsize * sizeof(p->t[0]));

    pthread_mutex_init(&p->pool_m, NULL);
    pthread_cond_init(&p->empty_c, NULL);
    pthread_cond_init(&p->pending_c, NULL);
    pthread_cond_init(&p->full_c, NULL);

    for (i = 0; i < tsize; i++) {
	if (0 != pthread_create(&p->t[i], NULL, t_pool_worker, p))
	    return NULL;
    }
 
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
    } else {
	j->serial = 0;
    }

#ifdef DEBUG
    fprintf(stderr, "Dispatching job %p for queue %p, serial %d\n", j, q, j->serial);
#endif

    pthread_mutex_lock(&p->pool_m);

    // Check if queue is full
    while (p->njobs == p->qsize)
	pthread_cond_wait(&p->full_c, &p->pool_m);

    p->njobs++;

    if (p->tail) {
	p->tail->next = j;
	p->tail = j;
    } else {
	p->head = p->tail = j;
    }

    if (p->njobs == 1) {
	// First job => tell all worker threads to start up
	pthread_cond_broadcast(&p->pending_c);
    }

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
	j->serial = q->curr_serial;
	pthread_mutex_unlock(&q->result_m);
    } else {
	j->serial = 0;
    }

#ifdef DEBUG
    fprintf(stderr, "Dispatching job for queue %p, serial %d\n", q, j->serial);
#endif

    pthread_mutex_lock(&p->pool_m);

    if (p->njobs == p->qsize && nonblock == 1) {
	pthread_mutex_unlock(&p->pool_m);
	errno = EAGAIN;
	free(j);
	return -1;
    }

    if (q) {
	pthread_mutex_lock(&q->result_m);
	q->curr_serial++;
	q->pending++;
	pthread_mutex_unlock(&q->result_m);
    }

    // Check if queue is full
    if (nonblock == 0)
	while (p->njobs == p->qsize)
	    pthread_cond_wait(&p->full_c, &p->pool_m);

    p->njobs++;
    
//    if (q->curr_serial % 100 == 0)
//	fprintf(stderr, "p->njobs = %d    p->qsize = %d\n", p->njobs, p->qsize);

    if (p->tail) {
	p->tail->next = j;
	p->tail = j;
    } else {
	p->head = p->tail = j;
    }

#ifdef DEBUG
    fprintf(stderr, "Dispatched (serial %d)\n", j->serial);
#endif

    if (p->njobs == 1) {
	// First job => tell all worker threads to start up
	pthread_cond_broadcast(&p->pending_c);
    }

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
#ifdef DEBUG
    fprintf(stderr, "Flushing pool %p\n", p);
#endif

    // Drains the queue
    pthread_mutex_lock(&p->pool_m);
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

	pthread_cond_broadcast(&p->pending_c);
	pthread_mutex_unlock(&p->pool_m);

#ifdef DEBUG
	fprintf(stderr, "Shutdown complete\n");
#endif
	for (i = 0; i < p->tsize; i++)
	    pthread_join(p->t[i], NULL);
    } else {
	for (i = 0; i < p->tsize; i++)
	    pthread_kill(p->t[i], SIGINT);
    }

    pthread_mutex_destroy(&p->pool_m);
    pthread_cond_destroy(&p->empty_c);
    pthread_cond_destroy(&p->pending_c);
    pthread_cond_destroy(&p->full_c);

#ifdef DEBUG_TIME
    fprintf(stderr, "Total time=%f\n", p->total_time / 1000000.0);
    fprintf(stderr, "Wait  time=%f\n", p->wait_time  / 1000000.0);
    fprintf(stderr, "%d%% utilisation\n",
	    (int)(100 - ((100.0 * p->wait_time) / p->total_time + 0.5)));
#endif

    free(p->t);
    free(p);

#ifdef DEBUG
    fprintf(stderr, "Destroyed pool %p\n", p);
#endif
}


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
    } else {
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
    t_pool *p = t_pool_init(NTHREADS*2, NTHREADS);
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
