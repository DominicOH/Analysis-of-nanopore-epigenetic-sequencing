from functools import wraps
import time

def timer(func):
    @wraps(func)
    def inner_func(*args, **kwargs):
        start = time.perf_counter()
        func_out = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"Elapsed time: {round(end-start, 3)} seconds.")
        return func_out
    return inner_func