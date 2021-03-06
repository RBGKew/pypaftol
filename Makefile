GITBRANCH       := $(shell git rev-parse --abbrev-ref HEAD)
GITTAG          := $(shell git describe --tags)


build : ptversion
	python setup.py build

hinstall : tests
	python setup.py install --home=$${HOME}

tests : build
	python setup.py test

ptversion :
	echo "branch: $(GITBRANCH), tag $(GITTAG)"
	echo "__version__ = '$(GITBRANCH)-$(GITTAG)'" > paftol/version.py

doc :
	epydoc -v paftol

clean :
	rm -rf html build dist paftol.egg-info paftol/version.py

.PHONY : build ptversion doc clean hinstall tests gitversion ptversion
