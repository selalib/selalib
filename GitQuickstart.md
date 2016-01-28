*Quickstart guide to use git and selalib, it's an email i sent to Morgane
with some additions. So i keep her name in it.* [Pierre](@pin)



Morgane wants to contribute to selalib project, she needs to

* Have an account on https://gitlab.mpcdf.mpg.de
* Send an email to @pin to be part of the project selalib
* After Pierre accept her in the projet, she has to edit her ssh keys
`https://gitlab.mpcdf.mpg.de/profile/keys`
and follow instructions...

After a couple of hours, she can install git and configure it
~~~
 git config --global user.name "Morgane Name"
 git config --global user.email "morgane@somewhere.fr"
 git config --global core.editor emacs (or other)
 git config --global merge.tool opendiff (or other)
~~~

Check configuration with:
```
 git config --list
```

Clone the repository:
~~~
 git clone git@gitlab.mpcdf.mpg.de:selalib/selalib.git
 cd selalib/
~~~

Display all branches with:
~~~
 git branch -a
~~~

Switch to the local default branch named develop:
~~~
 git checkout develop
~~~

To build the library, she reads the [CMake Quickstart Guide](CMakeQuickstart.md)
in this directory...

* To contribute it's better to create her own new branch
~~~
 git branch morgane-branch
~~~

* To switch to her new branch
~~~
 git checkout morgane-branch
~~~

* She wants to get all new updates of develop (the default branch)
~~~
 git fetch origin              
 git merge origin/develop
~~~

* She has conflict, no problem just do :
~~~
 git mergetool
~~~
It will open a nice editor and help her to choose the right version. Close the 
editor and commit.
~~~
 git commit -m 'Merge with selalib-devel and fixed conflicts'
~~~

* She implemented a new algorithm. To save it :
~~~
 git add new_algo_functions.F90
~~~

* She checks what files are ready to be committed, if there are untracked files, 
or if some tracked file has been modified but has not been staged for the next commit.
~~~
 git status
~~~

* She can now save her work on the local branch.
~~~
 git commit -m 'A new algorithm is available'
~~~

* If morgane-branch never pushed, it exists only on her disk.
But she wants to offer her new algorithm to the group.
She goes back to the main branch:
~~~
 git checkout develop   (switch to develop)
 git merge morgane-branch       (add the new algorithm from Morgane)
 git push origin develop  (develop is updated on the server, 
                                 everybody can use these improvements)
~~~

* But Morgane has an old computer and she wants to use the big hpc server with 
many cores, so she pushes her branch to the server
~~~
 git push origin morgane-branch
~~~

She logs on the cluster login node. She edit new ssh keys and configure git. To 
put her work on the server:
~~~
 git clone git@gitlab.mpcdf.mpg.de:selalib/selalib.git
 cd selalib/
 git checkout -t origin/morgane-banch
~~~