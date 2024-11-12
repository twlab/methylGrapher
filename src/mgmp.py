

import os
import sys
import time
import multiprocessing




# write a framework for both single threaded and multi-threaded execution

# 3 main class
# 1. Functions to execute
# 2. Worker class
# 3. Task class




# The single threaded execution is just one function after another
# The multi-processed execution is more complex



def test_step1(iter):

    i = 0
    while i < 10000:
        i += 1
        yield i


def test_step2(iter):

    for i in iter:
        time.sleep(0.0001)
        yield i**2





def worker_function(pid, name, function, input_queue, output_queue, batch_size=1024):
    print(f"{name} worker {pid} started")
    while True:
        try:
            tasks = input_queue.get()
            if tasks is None:
                break

            # print(pid, "tasks")

            result = []
            for r in function(tasks):
                result.append(r)

                if len(result) >= batch_size:
                    output_queue.put(result)
                    result = []

            if len(result) > 0:
                output_queue.put(result)

        except Exception as e:
            print(e)
            break

    print(f"{name} worker {pid} finished")


def end_function(pid, last_queue):
    print(f"End function {pid} started")
    while True:
        try:
            tasks = last_queue.get()
            if tasks is None:
                break

            # print(pid, len(tasks))
        except Exception as e:
            print(e)
            break
    print(f"End function {pid} finished")



class TaskOLD:

    def __init__(self):
        pass

    def run(self, threads=1):
        maxsize = 100
        batch_size = 256

        q0 = multiprocessing.Queue(maxsize=maxsize)
        q1 = multiprocessing.Queue(maxsize=maxsize)
        q2 = multiprocessing.Queue(maxsize=maxsize)



        pool1 = []
        for i in range(1):
            p = multiprocessing.Process(
                name=f"MGWorke1-{i}",
                target=worker_function,
                args=(i, "1", test_step1, q0, q1),
                kwargs={"batch_size": batch_size}
            )
            p.start()
            pool1.append(p)

        pool2 = []
        for i in range(threads):
            p = multiprocessing.Process(
                name=f"MGWorker2-{i}",
                target=worker_function,
                args=(i, "2", test_step2, q1, q2),
                kwargs={"batch_size": batch_size}
            )
            p.start()
            pool2.append(p)

        p = multiprocessing.Process(
            name=f"MGEnd",
            target=end_function,
            args=(0, q2)
        )
        p.start()




        q0.put(1)
        q0.put(None)

        for p in pool1:
            p.join()

        for p in pool2:
            q1.put(None)

        for p in pool2:
            p.join()

        q2.put(None)



        print("Finished")

        return None




# TODO how to pass parameters
class Task:

    def __init__(self, steps, batch_size=1024, queue_max_size=100):

        # Example steps
        steps_example = [
            {"name": "step1", "function": test_step1, "args": [], "kwargs": {}, "threads": 1},
            {"name": "step2", "function": test_step2, "args": [], "kwargs": {}, "threads": 2},
        ]

        self._steps = steps
        self._batch_size = batch_size
        self._queue_max_size = queue_max_size

        return

    def run(self):
        threads = 10
        maxsize = 100
        batch_size = 256


        # N steps need N queues, but for formatting purposes, it need N+1 queues
        q0 = multiprocessing.Queue(maxsize=self._queue_max_size)
        q1 = multiprocessing.Queue(maxsize=self._queue_max_size)
        q2 = multiprocessing.Queue(maxsize=self._queue_max_size)


        # N steps need N pools
        pool1 = []
        for i in range(1):
            p = multiprocessing.Process(
                name=f"MGWorke1-{i}",
                target=worker_function,
                args=(i, "1", test_step1, q0, q1),
                kwargs={"batch_size": batch_size}
            )
            p.start()
            pool1.append(p)

        pool2 = []
        for i in range(3):
            p = multiprocessing.Process(
                name=f"MGWorker2-{i}",
                target=worker_function,
                args=(i, "2", test_step2, q1, q2),
                kwargs={"batch_size": batch_size}
            )
            p.start()
            pool2.append(p)

        p = multiprocessing.Process(
            name=f"MGEnd",
            target=end_function,
            args=(0, q2)
        )
        p.start()





        # Let the steps begin
        q0.put([1])
        q0.put(None)

        for p in pool1:
            p.join()

        for p in pool2:
            q1.put(None)

        for p in pool2:
            p.join()

        q2.put(None)



        print("Finished")

        return None



def myfunc(x, y, z):
    print(x)
    print(y)
    print(z)
    print(x,y,z)
    return x + y + z


def test_args(a, b, c, *args, kw=3, **kwargs):
    print(a)
    print(b)
    print(c)
    print(kw)
    print(args)
    print(kwargs)
    myfunc(*args)
    return None




def to_worker_function(function, pid, name, input_queue, output_queue, *args, **kwargs):


    def wrapper(function, pid, name, input_queue, output_queue, *args, **kwargs):
        batch_size = kwargs.get("batch_size", 1024)
        print(f"{name} worker {pid} started")
        while True:
            try:
                tasks = input_queue.get()
                if tasks is None:
                    break

                # print(pid, "tasks")

                result = []
                for r in function(tasks):
                    result.append(r)

                    if len(result) >= batch_size:
                        output_queue.put(result)
                        result = []

                if len(result) > 0:
                    output_queue.put(result)

            except Exception as e:
                print(e)
                break

        print(f"{name} worker {pid} finished")

    return wrapper





test_args(1, 2, 3, 100, 200, 300, verbose=True, silent=False)


sys.exit(3)

if __name__ == "__main__":


    # Single threaded execution

    ts1 = time.time()
    s1_gen = test_step1(lambda x: x)
    s2_gen = test_step2(s1_gen)
    for i in s2_gen:
        pass
    duration1 = time.time() - ts1


    # Multi-threaded execution

    ts2 = time.time()
    task = Task(1)
    task.run()
    duration2 = time.time() - ts2

    ts5 = time.time()
    task = Task(1)
    task.run()
    duration5 = time.time() - ts5

    ts10 = time.time()
    task = Task(1)
    task.run()
    duration10 = time.time() - ts10

    print(duration1, duration2, duration5, duration10)
























