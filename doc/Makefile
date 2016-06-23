default: doc html

html:
	(cd sphinx; make html)
sync: html
	rsync --delete -vurz -e ssh ./build/html/* scm.gforge.inria.fr:/home/groups/selalib/htdocs/ 
doc:
	(cd ../build; make doc-user)
	#(cd ../build; make doc-dev)

