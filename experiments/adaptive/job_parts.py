"""Standard Job worker parts with values retained at stage completion."""

from contextlib import contextmanager
import time

from microgen_workflow import MEMORY


@contextmanager
def worker_part(job, name):
    part = job.part(name)
    started = time.monotonic()
    part.update(0, "starting bounded worker", force=True, elapsed_s=0, memory_ceiling_gib=MEMORY / 1024**3)
    try:
        yield part
    except BaseException as error:
        part.failed(f"failed: {type(error).__name__}", elapsed_s=time.monotonic() - started,
            memory_ceiling_gib=MEMORY / 1024**3)
        raise
    else:
        part.done("done", elapsed_s=time.monotonic() - started, memory_ceiling_gib=MEMORY / 1024**3)
