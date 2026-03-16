import icikt
import time
import numpy

in_data = numpy.random.randn(10000, 400)

start_time = time.perf_counter()

# Code to be timed
scaled, corrRaw, pVals, tauMax = icikt.iciktArray(in_data)
# for i in range(10000000):
#     pass

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.4f} seconds")
