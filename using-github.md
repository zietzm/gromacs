# Using GitHub
We are using GitHub for this project so that we can use native text editors on
Windows and Ubuntu. However, this process is not as seamless as we would hope.
<br><br>
Whenever changes are made at GitHub or locally on any machine, we have to be
<b>very</b> careful that we do not get ahead of or behind other versions of
ourselves. This is definitely doable, but it will require some discipline and
care.
<br><br>
To begin on a new machine, we would need to clone our Git repository from
GitHub.
<pre><code>git clone git@github.com:zietzm/gromacs-information.git
 </code></pre>

Now, an outline of what we have to do when we are working locally.
1. Save our altered file locally. Ctrl-s works for this.
2. Stage our altered file.
<pre><code>git add FILE.EXT</pre></code>
3. Commit changes.
<pre><code>git commit -m "Comment" </pre></code>
4. Push changes to GitHub.
<pre><code>git push origin master</pre></code>
