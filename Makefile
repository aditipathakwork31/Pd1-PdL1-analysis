# Single documented run sequence. Each target is added as the step is completed.
PY ?= .venv/bin/python

env:
	python3 -m venv .venv
	.venv/bin/pip install -r requirements.txt
	.venv/bin/pip freeze > requirements.lock
	$(PY) -c "import sys,Bio,numpy,platform;print(platform.python_version(),Bio.__version__,numpy.__version__)" | tee outputs/versions.txt

fetch:
	$(PY) src/fetch_inputs.py fetch

verify:
	$(PY) src/fetch_inputs.py verify

compare:
	$(PY) src/fetch_inputs.py compare --dir $(DIR)

link:
	$(PY) src/fetch_inputs.py link

test:
	$(PY) -m pytest -q tests

.PHONY: env fetch verify link compare test
