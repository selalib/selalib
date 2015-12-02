default: doc html

sync:
	rsync --delete -vurz -e ssh ./build/html/* scm.gforge.inria.fr:/home/groups/selalib/htdocs/ 
html:
	(cd sphinx; make html)
doc:
	(cd ../build; make doc-user)
	#(cd ../build; make doc-dev)

