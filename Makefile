default:
	python setup.py build

hinstall : tests
	python setup.py install --home=$${HOME}

tests :
	python setup.py test

doc :
	epydoc -v paftol

clean :
	rm -rf html build dist paftol.egg-info

.PHONY : default doc clean hinstall tests
