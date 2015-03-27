all: cythonize

cythonize: fast_metrics.pyx
	python setup.py build_ext --inplace

clean:
	rm -rf *c *o
