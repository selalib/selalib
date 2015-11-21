Quickstart guide to use git and selalib, it's an email i sent to Morgane
with some additions. So i keep her name in it.

Pierre

Morgane wants to contribute to selalib project, she needs to

* Create an account on gforge.inria.fr
* Send an email to navaro@math.unistra.fr to be part of the project selalib
* After Pierre accept her in the projet, she has to edit her ssh keys
`https://gforge.inria.fr/account/editsshkeys.php`
and follow instructions...

After a couple of hours, she can install git and configure it
~~~
 git config --global user.name "Morgane Name"
 git config --global user.email "morgane@somewhere.fr"
 git config --global core.editor emacs (or other)
 git config --global merge.tool opendiff (or other)
~~~

Check configuration with:
 git config --list

Clone the repository:
~~~
 git clone git+ssh://morgane@scm.gforge.inria.fr//gitroot//selalib/selalib.git
 cd selalib/
~~~

Display all branches with:
~~~
 git branch -a
~~~

Create the local branch and merge the remote selalib-devel branch:
~~~
 git checkout -b selalib-devel origin/selalib-devel
~~~

To build the library, she reads the CMakeQuickstart.txt in this directory...

* To contribute it's better to create her own new branch
~~~
 git branch morgane-branch
~~~

* To switch to her new branch
~~~
 git checkout morgane-branch
~~~

* She wants to get all new updates of selalib-devel (the main branch)
~~~
 git checkout selalib-devel      (change to selalib-devel)
 git fetch origin                  (download changes from repository)
 git merge origin/selalib-devel  (update local branch selalib-devel)
 git checkout morgane-branch       (back to her branch)
 git merge selalib-devel         (update her branch with modifications on selalib-devel)
~~~

* She has conflict, no problem just do :
~~~
 git mergetool
~~~
It will open a nice editor and help her to choose the right version. Close the editor and commit.
~~~
 git commit -m 'Merge with selalib-devel and fixed conflicts'
~~~

*She implemented a new algorithm. To save it :
~~~
 git add new_algo_functions.F90
~~~

*She checks what files are ready to be committed, if there are untracked files, or if some tracked file has been modified but has not been staged for the next commit.
~~~
 git status
~~~

*She can now save her work on the local branch.
~~~
 git commit -m 'A new algorithm is available'
~~~

*If morgane-branch never pushed, it exists only on her disk (invisible) ~~~~~(m--)m
But she wants to offer her new algorithm to the group.
She goes back to the main branch:
~~~
 git checkout selalib-devel   (switch to selalib-devel)
 git merge morgane-branch       (add the new algorithm from Morgane)
 git push origin selalib-devel  (selalib-devel is updated on the server, everybody can use these improvements)
~~~

*But Morgane has an old computer and she wants to use the big hpc server with many cores, so she pushes her branch to the server
~~~
 git push origin morgane-branch
~~~

She logs on the cluster login node. She edit new ssh keys and configure git. To put her work on the server:
~~~
 git clone git+ssh://morgane@scm.gforge.inria.fr//gitroot//selalib/selalib.git
 cd selalib/
 git checkout -b morgane-branch origin/morgane-banch
~~~