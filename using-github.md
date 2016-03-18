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

<b>Now, an outline of saving local work to GitHub.</b><br><br>
1. Save our altered file locally. Ctrl-s works for this.<br>
2. Stage our altered file. <pre><code>git add FILE.EXT</code></pre>
3. Commit changes. <pre><code>git commit -m "Comment" </code></pre>
4. Push changes to GitHub. <pre><code>git push origin master</code></pre>
