You want to contribute to selalib project, you need to

* Have an account on https://gitlab.mpcdf.mpg.de
* Edit your ssh keys
`https://gitlab.mpcdf.mpg.de/profile/keys`
and follow instructions...

After a couple of hours, install git and configure it
~~~
 git config --global user.name "Your Name"
 git config --global user.email "your.name@somewhere.fr"
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

To build the library, read the [CMake Quickstart Guide](CMakeQuickstart.md)

* To contribute it's better to create your own branch
~~~
 git branch your-branch
~~~

* Switch to the new branch
~~~
 git checkout your-branch
~~~

* To get all new updates of develop (the default branch)
~~~
 git fetch origin              
 git rebase origin/develop
~~~

* If you have conflicts, no problem just do :
~~~
 git mergetool
~~~
It will open a nice editor and help you to choose the right version. Close the 
editor and commit.
~~~
 git commit -m 'Merge with selalib-devel and fixed conflicts'
~~~

* You implemented a new algorithm. To save it :
~~~
 git add new_algo_functions.F90
~~~

* Check what files are ready to be committed, if there are untracked files, 
or if some tracked file has been modified but has not been staged for the next commit.
~~~
 git status
~~~

* Save your work on the local branch.
~~~
 git commit -m 'A new algorithm is available'
~~~

* If your-branch never pushed, it exists only on your disk.
To offer the new algorithm to the group. Push your branch on the server:
~~~
 git push origin your-branch
~~~
And (file a merge request)[https://gitlab.mpcdf.mpg.de/clapp/selalib/merge_requests].

If you want to change of computer, you can get your branch with:
~~~
 git clone git@gitlab.mpcdf.mpg.de:selalib/selalib.git
 cd selalib/
 git checkout -t origin/your-branch
~~~

When the work is done, you can delete the remote branch:
~~~
git push origin --delete your-branch
~~~

You can also remove the local branch by switching to another branch and delete it:
~~~
git checkout master
git branch -D your-branch
~~~
