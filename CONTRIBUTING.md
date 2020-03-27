# Introduction

NEPTUNE is used in many IRAS software projects. Developers of such software 
projects normally want to use and modify NEPTUNE. They add functionality to it, 
improve existing functions and fix bugs. The problem:

1. They don't want to be bothered by other developers using NEPTUNE.
2. They want to use all the updates other developers using NEPTUNE do.

This description shall bring order into this mess. It is based on the process
setup by SvM for libslam!

# Branches

On github, there is one **master branch** named "master" which contains the most 
recent version of NEPTUNE. Each separate change to NEPTUNE shall be done in its 
own **change branch** named "<Developer Abbreviation>/<Change Name>" where 
"<Devloper Abbreviation> must be your official abbreviation in lower-case 
(e. g. "jr" or "anh") and "<Change Name>" must be something like "fix_nxtbuf" 
or "add_nasabum". Please always use underscores ("_") for word separation.

Each branch has at least one owner. Owner of the master branch is ChK. Owner of 
a change branch is the developer working on the change which is supposed to be 
introduced. The developer abbreviation in a change branch name indicates the 
owner of the branch. **Only the owner of a branch may push changes to his/her 
branch on gitlab.**

# Making a Change

## Working in your Change Branch

You can make a change in the world (of NEPTUNE)! Here is what you have to do 
- assuming you are "ChK" and you want to "fix_nxtbuf":

Clone neptune git:

```
mkdir neptune-git
cd neptune-git
git clone http://github.com/Space-Systems/neptune.git .
```


If you use an already-existing local git repository instead:

```
cd neptune-git
git checkout master
git pull --ff-only
```

Note that you are, in both cases, in the master branch now. Create a new change 
branch based on master:

```
git branch chk/fix_nxtbuf
git checkout chk/fix_nxtbuf
```

NOTE: You are encouraged to create change branches for each change that you want 
to introduce to NEPTUNE - even if it is a small one. The smaller the change, the 
easier the merge at the end. As a result: The smaller the change, the sooner 
your change will arrive in NEPTUNE master.

If you are working on an "official" NEPTUNE issue in gitlab, you can use the 
"New branch" button there to create your change branch off the issue. However, 
the button will directly create a branch with the issue's number and title in 
its name (You can't choose another name). The title of an issue does not always 
name the change that actually needs to be done and there could be more than one 
change necessary to close the issue. Also, your abbreviation is missing in the 
branch name. So, this is bullocks. Interestingly though, gitlab will ask for the 
name of the new branch if a branch with the automatically-generated name already 
exists. So, you are welcome to:
1. Click on "New branch"
2. Go back to the issue
3. Click on "New branch" again
4. Change the branch name
5. Go to the branch with the "wrong" name
6. Delete the branch with the "wrong" name

After this you can clone the NEPTUNE git as shown above and then checkout your 
branch:

```
git checkout svm/3_fix_nxtbuf
```

Or, if you are using an already-existing repository:

```
git fetch
git checkout svm/3_fix_nxtbuf
```

Supposedly, you already noted the "3" in the branch name. This is the issue 
number of this example. Please always put your issue number into the name as 
shown above (if you are working on an "official" issue).

Now, let's continue. Make your changes and commit them:

```
(Change, for example, slam_io.f90)
git add slam_io.f90
git commit
```
**All commit messages must be prefixed by "svm/fix_nxtbuf: "**

**Always be compliant to NEPTUNE's [code conventions](Code conventions)**. 
If you are not, you will have to make your code compliant later anyhow. So, 
it's always better to be compliant from the start. The more compliant you are 
from the beginning, the sooner your changes will be merged into master.

NOTE: Each commit should cover a (small) change. Please don't put two different 
changes into one and the same commit. The commit message should summarize the 
change.

Afterwards, push your branch to github:

```
git push --set-upstream origin chk/fix_nxtbuf
```

After your inital push, you can do further pushes with only:

```
git push
```

NOTE: If you pushed before, but it doesn't work anymore it may be the case that 
your branch has already been merged into master and then removed from github! 
In this case, if something is missing, please follow the instructions under 
"Quick-Fixing" beneath.

You can make as many adds, commits and pushes as you like and over a period of 
time as long as you like. You can also work on more than one change in parallel. 
In any case, it is your responsibility to always stay up-to-date with master! 
To accomplish this, you need to **rebase on master**. This means the following:

1. Your commits on master are temporarily put aside.
1. Your change branch is set equal to master.
1. Each of your commits is transformed into a new commit on master one after 
another. This new commit will be the same as the old one by default. But in the 
case of conflicts, you will need to resolve the conflict which makes the new 
commit different from the old one.

In simple words, what you want to do by rebasing is to introduce your changes 
one by one again, based on the new version of master. So, you rebase your 
branch on the (new) top of master. This will lead to a cleaner history in the 
end. Translated into commands, this means you have to first do the rebase 
itself:

```
git checkout master
git pull --ff-only
git checkout chk/fix_nxtbuf
git rebase master
```

If a conflict occurs, you will be asked to resolve it. In the conflicted file 
you will find something like:

```
configure_file(src/slam.h.in src/slam.h) # COMMENT FROM CHANGE BRANCH chk/fix_nxtbuf
```

In this case a comment beginning with "#" was added in another change branch 
named "svm/something" to the same code line in CMakeLists.txt. Note that you 
will normally not see the branch name of the other change branch. It has been 
added here to the comment for explanatory purposes. "HEAD" is the new HEAD of 
your change branch = master HEAD + commits put ontop of it up to the current 
commit which produced the conflict. You now need to resolve the conflict like 
normal and then continue with the rebase. For example:

```
(Change CMakeLists.txt)
git add CMakeLists.txt
git rebase --continue
```

Here, it is possible that git tells you that you can't continue, since no 
changes have been commited. In this case, one possible explanation is that you 
erased your own changes made to CMakeLists.txt in the currently-processed 
commit. If that is the case and it was intended as such, use the following 
command instead of the continue command to proceed:

```
(Change CMakeLists.txt)
git add CMakeLists.txt
git rebase --skip
```

NOTE: If, after a rebase, you are adviced by git to pull changes: Don't. It's a 
trap!

After the rebase, you need to push your changes to github. A simple push will not 
do the trick this time, but the following will:

```
git push origin +chk/fix_nxtbuf
```

Note the "+"! This makes this push a force push which overwrites your branch 
on github completely - including its history. DO NOT use "--force" or "-f" 
instead, since this may lead to a force push of all branches!

## Creating a Merge Request


As soon as you think your change is ready to be incorporated into master, you 
need to make a **merge request**:

1. Check if your code is compilable and conforms to the [code conventions]
(Code conventions).
1. Open github's github web interface under http://github.com/Space-Systems and login.
1. Click on "neptune" project.
1. Click on "branches".
1. Click on "Merge request" next to your change branch.
1. Set "Asignee" to master owner.
1. Click on "Submit merge request".

With this, the master owner should be notified by a brown marker at the 
top-right of github's web interface. He/She will:

1. Look at your commits.
1. Pull and checkout your change branch
1. Compile
1. Run the tests.

## Merge Request Acceptance

If everything is fine, the master owner will merge your change branch into 
master as-is. In this case:

1. The master owner will notify you about the merge.
1. The master owner creates a new tag, thereby increasing NEPTUNE's version 
number.
1. The master owner removes your change branch from github.
1. You must remove your change branch from your local system/repository.

## Merge Request Iteration

If, however, your change needs... changes or clarification, the master owner 
will:

1. Add comments to your merge request.
1. Set the assignee of the merge request to you. This will make the brown merge 
request notifier pop up in your own gitlab view.

Then, please:

1. Open the merge request yourself in gitlab.
1. You will be able to see the comments by scrolling down there. Read the 
comments.

You are able to comment on your merge request yourself. If you think, the next 
step should be for the master owner to look at your comments:

1. In the merge request, click on "Edit".
1. Change the "Assignee" back to the master owner.

If now, your merge request is accepted by the master owner, see above. In any 
case, you can (and - depending on the comments - should):

1. Make further adds, commits and pushes (and rebases!) to your change branch. 
These new changes will be added to your existing merge request. You don't have 
to create a new merge request.
1. Set the "Assignee" of the merge request (back) to the master owner.

NOTE: It is always possible that you will be asked to do rebases on master by 
the master owner, before your request can be processed! This will happen as 
soon as a merge request is accepted before your's.

## Quick-Fixing

If you become aware that your merge request has been accepted, but you forgot 
something OR if you just want to fix something quickly, you have to do it in a 
new change branch following the same process as above. However, it doesn't take 
as long as you might think:

```
git checkout master
git pull --ff-only
git branch chk/quickfix_sth
git checkout chk/quickfix_sth
git add ...
git commit
git push --set-upstream origin chk/quickfix-sth
```

Then, create a merge request for your change branch and set the assignee to the 
master owner. Done.

TODO: Does commiting to an already-merged branch, rebasing it and pushing work, 
too?

## Temporary Branches

Your merge request will be processed as soon as possible. If you or your project 
needs a new version of NEPTUNE sooner as the necessary merge requests are 
processed, you should create a **temporary branch**. If, for example, you need a 
NEPTUNE version for RSS Project Meeting 8:

```
git checkout master
git branch tmp_neptune_rss8
git checkout chk/fix_nxtbuf
git pull
git checkout sh/do_something
git pull
git checkout tmp_neptune_rss8
git merge chk/fix_nxtbuf
git merge sh/do_something
git push --set-upstream origin tmp_neptune_rss8
```

If possible, you should always keep working in your own change branches, though, 
and not in your temporary branch:

```
git checkout chk/fix_nxtbuf
git add ...
git commit
git push
git checkout tmp_neptune_rss8
git merge chk/fix_nxtbuf
git push
```