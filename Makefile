default:
	python setup.py build

hinstall :
	python setup.py install --home=$${HOME}

doc :
	epydoc paftol

clean :
	rm -rf html build dist paftol.egg-info

.PHONY : default doc clean hinstall

