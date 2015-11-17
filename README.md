Dear selalib developers,


To work with gitlab you don’t need to clone the repository and stay with
the clone you download from INRIA gorge. If you want to push a branch on 
gitlab repository, in your selalib directory add the repository with
$ git remote add gitlab git@gitlab.mpcdf.mpg.de/selalib/selalib.git

After this, you can push your local branch on gitlab with
$ git push gitlab your_branch

Your branch will be available on gitlab interface and available for a merge request
If this branch is changed on the interface, you can update it with:
```bash
    git fetch gitlab
    git checkout your_branch
    git merge gitlab/your_branch
```

You can push to the INRIA gorge
$ git push origin your_branch

If you merge your branch with the develop branch on gitlab, don’t forget to remove
your_branch on gitlab.
$ git push gitlab —delete your_branch

it deletes the remote branch not the local one and you still have a copy on the INRIA gforge.

With this new gitlab server, we change the way we use git.
I describe this process, you can find the long version with figures here
 (https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow)
I added some git commands but branch creation and merge could be done through the gitlab interface.


New rules

- Developers don't freely merge into main branch.
- They must file a merge request asking to merge their additions. 
- The new workflow assigns very specific roles to different branches.

New Branches Management

- We will uses two branches to record the history of the project. 
	- The master branch stores the official release history.
	- The develop branch serves as an integration branch for features.

Every new branch created has to be dedicated to a new feature that will be 
available in the release. If not, you can stay with the actual repository on the INRIA gforge.
The master branch on the INRIA gforge will be regularly updated from master branch of the gitlab server.

- Each new feature should reside in its own branch. 
- Feature branches use develop as their parent branch. 
- Features should never interact directly with master.

Example: Develop a new feature

Base the feature branches on develop:
$ git checkout -b some-feature develop
After some commits, file a merge request and delete the branch:
$ git pull gitlab develop
$ git checkout develop
$ git merge some-feature
$ git push gitlab develop
$ git branch -d some-feature

Very important, the branch created has to be deleted after the merge. If you want to
fix something or extend a feature, create a new branch!

Release Process

- Project owner creates a release branch off of develop. 
- No new features can be added after this point, only bug fixes, documentation. 
- Once it's ready to ship, the release gets merged into master. 
- In addition, it should be merged back into develop, (which may have progressed).
- The release process concerns everybody because documentation and bug fix
cannot be done only by project owners. 

The following commands will be done only by one person. 
Every achieved step has to be communicated to all developers.

Example: Prepare the release

$ git checkout -b release-0.1 develop
After documentation generation and bug fixes, merge it into master and develop.
$ git checkout master
$ git merge release-0.1
$ git push gitlab master
$ git checkout develop
$ git merge release-0.1
$ git push gitlab develop
$ git branch -d release-0.1

Tag the commit :
$ git tag -a 0.1 -m "Initial public release" master
$ git push --tags

For Bug Fixing

This part concerns every developer and bug should be reported and assigned through the “Issues” tab of Gitlab.
Create a maintenance branch off of master, fixes the issue with as many commits as necessary, 
then merges it directly back into master.

$ git checkout -b issue-#001 master
  Fix the bug
$ git checkout master
$ git merge issue-#001
$ git push gitlab master

Updates  need to be included also in develop

$ git checkout develop
$ git merge issue-#001
$ git push gitlab develop
$ git branch -d issue-#001
