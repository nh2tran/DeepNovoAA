clean:
	rm -rf build
	rm -f deepnovo_cython_modules.c
	rm -f deepnovo_cython_modules*.so

pull:
	git pull

.PHONY: build
build: clean
	python deepnovo_cython_setup.py build_ext --inplace

.PHONY: train
train: pull
	python deepnovo_main.py --train

.PHONY: denovo
denovo:
	python deepnovo_main.py --search_denovo

.PHONY: denovo_test
denovo_test: train denovo
	python deepnovo_main.py --test
