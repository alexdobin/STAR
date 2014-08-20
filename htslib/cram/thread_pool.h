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

/*
 * This file implements a thread pool for multi-threading applications.
 * It consists of two distinct interfaces: thread pools an results queues.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling.
 *
 * Upon completion, the return value from the function pointer is added to
 * a results queue. We may have multiple queues in use for the one pool.
 *
 * An example: reading from BAM and writing to CRAM with 10 threads. We'll
 * have a pool of 10 threads and two results queues holding decoded BAM blocks
 * and encoded CRAM blocks respectively.
 */

#ifndef _THREAD_POOL_H_
#define _THREAD_POOL_H_

#include <pthread.h>

struct t_pool;
struct t_results_queue;

typedef struct t_pool_job {
    void *(*func)(void *arg);
    void *arg;
    struct t_pool_job *next;

    struct t_pool *p;
    struct t_results_queue *q;
    int serial;
} t_pool_job;

typedef struct t_res {
    struct t_res *next;
    int serial; // sequential number for ordering
    void *data; // result itself
} t_pool_result;

typedef struct t_pool {
    int qsize;    // size of queue
    int njobs;    // pending job count
    int nwaiting; // how many workers waiting for new jobs
    int shutdown; // true if pool is being destroyed

    // queue of pending jobs
    t_pool_job *head, *tail;

    // threads
    int tsize;    // maximum number of jobs
    pthread_t *t;

    // Mutexes
    pthread_mutex_t pool_m; // used when updating head/tail

    pthread_cond_t  empty_c;
    pthread_cond_t  pending_c; // not empty
    pthread_cond_t  full_c;

    // Debugging to check wait time
    long long total_time, wait_time;
} t_pool;

typedef struct t_results_queue {
    t_pool_result *result_head;
    t_pool_result *result_tail;
    int next_serial;
    int curr_serial;
    int queue_len;  // number of items in queue
    int pending;    // number of pending items (in progress or in pool list)
    pthread_mutex_t result_m;
    pthread_cond_t result_avail_c;
} t_results_queue;


/*
 * Creates a worker pool of length qsize with tsize worker threads.
 *
 * Returns pool pointer on success;
 *         NULL on failure
 */
t_pool *t_pool_init(int qsize, int tsize);

/*
 * Adds an item to the work pool.
 *
 * FIXME: Maybe return 1,0,-1 and distinguish between job dispathed vs
 * result returned. Ie rather than blocking on full queue we're permitted
 * to return early on "result available" event too.
 * Caller would then have a while loop around t_pool_dispatch.
 * Or, return -1 and set errno to E_AGAIN to indicate job not yet submitted.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int t_pool_dispatch(t_pool *p, t_results_queue *q,
		    void *(*func)(void *arg), void *arg);
int t_pool_dispatch2(t_pool *p, t_results_queue *q,
		     void *(*func)(void *arg), void *arg, int nonblock);

/*
 * Flushes the pool, but doesn't exit. This simply drains the queue and
 * ensures all worker threads have finished their current task.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int t_pool_flush(t_pool *p);

/*
 * Destroys a thread pool. If 'kill' is true the threads are terminated now,
 * otherwise they are joined into the main thread so they will finish their
 * current work load.
 *
 * Use t_pool_destroy(p,0) after a t_pool_flush(p) on a normal shutdown or
 * t_pool_destroy(p,1) to quickly exit after a fatal error.
 */
void t_pool_destroy(t_pool *p, int kill);

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
t_pool_result *t_pool_next_result(t_results_queue *q);
t_pool_result *t_pool_next_result_wait(t_results_queue *q);

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
void t_pool_delete_result(t_pool_result *r, int free_data);

/*
 * Initialises a results queue.
 *
 * Results queue pointer on success;
 *         NULL on failure
 */
t_results_queue *t_results_queue_init(void);

/* Deallocates memory for a results queue */
void t_results_queue_destroy(t_results_queue *q);

/*
 * Returns true if there are no items on the finished results queue and
 * also none still pending.
 */
int t_pool_results_queue_empty(t_results_queue *q);

/*
 * Returns the number of completed jobs on the results queue.
 */
int t_pool_results_queue_len(t_results_queue *q);

/*
 * Returns the number of completed jobs plus the number queued up to run.
 */
int t_pool_results_queue_sz(t_results_queue *q);

#endif /* _THREAD_POOL_H_ */
