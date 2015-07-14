default: doc html
	rsync -vurz --exclude 'releases' --delete -e ssh ./build/html/* \
               scm.gforge.inria.fr:/home/groups/selalib/htdocs/ 
html:
	(cd sphinx; make html)
doc:
	(cd ../build; make doc-user)
	(cd ../build; make doc-dev)

