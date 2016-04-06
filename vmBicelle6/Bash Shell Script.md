# Working with Shell Scripts and Bash Shell
To make simulations easier and manual, we want to make them run automatically.
To this end, we can use shell scripts to run a series of commands without user
input.

We can write this with VIM in the terminal, NANO, or any other text editor we
prefer. If we are using VIM, for example, we type `vim FILENAME.sh`, where we
obviously replace "FILENAME" with the name of the file we wish to create.

Inside vim, press `i` to begin inserting (typing). Next, we need to type
```
#!/bin/bash
```
After this, we can type out any commands that run in the terminal. For just
commands, there is no need to even put semicolons at the end of text. If we
want to print output as well, however, we need to use `echo "TEXT";` to print
our desired "TEXT". Once we have finished writing all the commands that will
run automatically, simply press ESCAPE, followed by `:w` followed by `:x` to
save our work and close VIM. 
