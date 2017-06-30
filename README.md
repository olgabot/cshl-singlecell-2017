# cshl-singlecell-2017
Single Cell Analysis course at Cold Spring Harbor Laboratory 2017

This is one of *many* single cell courses/tutorials. An excellent list of all
single cell package, courses, tutorials, speakers for conferences, can be found
[here](https://github.com/seandavi/awesome-single-cell).

## Sources of Data

We will be using pre-cleaned data for the demos in the course. If you are
interested in how the data was cleaned (removing spaces, adding `"r1_"` to the
cell barcodes, combining multiple files, etc). Here are the GitHub repositories
for you to browse:

- [Macosko 2015](https://github.com/olgabot/macosko2015)
- [Segerstolpe + Palasantza et al 2016](https://github.com/olgabot/segerstolpe_palasantza2016)
- [Lönnberg + Svensson et al 2017](https://github.com/olgabot/lonnberg_svensson2017)


## Getting Started (email to students)

### Subject: Slack + Bioinformatics(Papers + software)

Dear Soon-to-be Single Cell Analysts,<br> I hope you are as excited as I am
about the upcoming class! We start in <2 weeks!! This is quite the beefy email
so please email us (the bioinformatics instructors: Olga, Emily, Alain) if you
have any questions. <br>
Warmest, <br>
Olga, Emily, Alain

### Slack messaging service
To help you get to know each other, I used the Slack messaging service to create a team for us: link (expires in one week - email me if you need access). You can sign in with your email and create a password. This is like a chat room that is only accessible to us! You can post pictures, papers, whatever you like. We use Slack at the CZ Biohub and for the Human Cell Atlas and it’s very helpful to stay in touch across the world! Before the course, I will be online as much as I can and will try to answer your questions as much as possible. During the course, it will be very useful to communicate with each other.


**Log on to Slack and introduce yourself in the #general channel by answering
the question: *“What song best represents your research and why?”* For example, I
picked Beyonce - “Single Ladies (Put a Ring on It)” because if you liked it,
you should have put a seq on it.**


### Papers

In the bioinformatics section, we will be performing case studies of four papers that embody the types of questions asked in single-cell RNA-seq analyses (and, by extension, single-cell genomics, epigenomics, proteomics, etc). To help you prepare, please read the following *papers*. If you are especially interested in the topic, there is also one or more *optional review(s)* that I think are excellent if you want more information.


During the course, we will be exploring their data so no need to download or analyze any data before the course.


1. **“Fishing expedition:” Take a tissue, dissociate it → ?? → subpopulations!**
   1. **Paper:** Macosko - 2015 - “Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets” (aka “the dropseq paper”)
https://www.ncbi.nlm.nih.gov/pubmed/26000488
   2. **Optional Reviews:**
Bacher and Kendziorski - 2016 - “Design and computational analysis of single-cell RNA-sequencing experiments”
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y
Kolodziejczyk et al - 2015 - “The Technology and Biology of Single-Cell RNA Sequencing”
http://www.cell.com/molecular-cell/fulltext/S1097-2765(15)00261-0
1. **Case vs Control**
   1. **Paper:** Segelstorpe - 2016 - “Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes” http://www.cell.com/cell-metabolism/fulltext/S1550-4131(16)30436-3
   2. **Optional Review:**
        Stegle et al - 2015 - “Computational and analytical challenges in single-cell transcriptomics”
        https://www.ncbi.nlm.nih.gov/pubmed/25628217
3. **Molecular transformations over time (“pseudo-time”)**
   1. **Paper:** Lönnberg et al -  2017  - “Single-cell RNA-seq and computational analysis using temporal mixture modelling resolves Th1/Tfh fate bifurcation in malaria”
https://www.ncbi.nlm.nih.gov/pubmed/28345074
   2. **Optional Reviews:**
Canoodt et al - 2016 - “Computational methods for trajectory inference from single-cell transcriptomics”
http://onlinelibrary.wiley.com/doi/10.1002/eji.201646347/epdf
Symmons and Raj - 2016 - “What’s Luck Got to Do with It: Single Cells, Multiple Fates, and Biological Nondeterminism”
https://www.ncbi.nlm.nih.gov/pubmed/27259209
4. **Perturbation**
    1. **Paper:**
    Adamson and Norman et al - 2016 - “A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response”
    https://www.ncbi.nlm.nih.gov/labs/articles/27984733/
    2. **Optional Review:**
    Tanay and Regev - 2017 - “Scaling single-cell genomics from phenomenology to mechanism” https://www.ncbi.nlm.nih.gov/pubmed/28102262
    Wagner and Klein - 2017 - “Genetic screening enters the single-cell era”
    http://www.nature.com/doifinder/10.1038/nmeth.4196


#### Reading questions

For each of these papers, answer the following
questions:


* What was their biological question?
* What were their quality control methods? Specifically ….
    * How did they filter cells?
    * How did they filter genes?
* How do they account for the following confounding factors (and any others)?
    * Cell size
    * Dead cells
    * Extremely high expression (e.g. insulin is very highly expressed in certain pancreas cell types)
    * Dropout of lowly expressed genes (e.g. transcription factors)
* For each of the following, answer what is the heterogeneity and how they accounted for it.
    * What is the known biological heterogeneity? (e.g. knockout of a gene)
    * What is the unknown biological heterogeneity? (e.g. the cell type heterogeneity)
    * What is the biological stochasticity? (e.g. transcriptional bursting)
    * What is the technical heterogeneity? (e.g. dropout or batch effects)
* What computational steps did they take to answer their biological question? It’s okay if these terms are over your head right now, we’ll go over them in the course. What’s important is to capture the overall flow chart.
    * E.g. Filtered genes → PCA → k-means clustering → differential expression

### Software
**Please make sure to do this before the course starts!**


We will be using primarily Python, because besides being the easiest language to teach, Python has a rich machine learning ecosystem that we will be taking advantage of, called scikit-learn.


However, due to the complexity of the low-level linear algebra libraries necessary for numerical Python, I do not recommend using the Python that you used in the Genomics course because you will likely have to do a lot of manual work to get them to work. Instead, we'll use the Anaconda Python Distribution because it comes pre-baked with a ton of helpful scientific packages for dataframes, matrices, statistics, machine learning, etc.


**If you already have Anaconda Python installed, jump to the section marked
between the <span style="background-color:cyan">blue highlighted START and
STOP</span>, if not, then no worries! Start right here at the section marked
between the <span style="background-color:yellow">yellow highlighted START and STOP</span>,and stop before you get to the
blue.**

<span style="background-color:yellow">**START here if you do not have Anaconda
Python installed.**</span>

1. **If you don't already have the Anaconda Python Distribution
(easiest!)**
   1. Please download and install the Anaconda Python Distribution for your laptop
   2. We need to install a few additional packages that don't already come with
      Anaconda, namely seaborn for plots. First you must locate your terminal
      program.
      1. (Windows) After installing, search for "Anaconda" in your windows
         search bar and click on "Anaconda command prompt".
      2. (Mac/Linux) Search for "Terminal" in the spotlight/program search and
         click on it.
    2. Once a terminal of some variety is open, type `conda install seaborn`
    1. Which will then prompt you for whether you're sure you want to install these.
   Type `y`(for yes) and press ENTER.

<span style="background-color:yellow">**STOP here if you started without Anaconda
Python - you do not need to make any environments!!**</span>


<span style="background-color:cyan">**START here if you already have Anaconda
Python installed.**</span>

**If you already have the Anaconda Python Distribution (harder)**

1. If it's already Python 3, you don't need to do anything (and can stop
   reading here), unless you want to be fancy and create a separate environment
   for this course. In that case, follow the instructions in 2.2.2 for creating
   a new Jupyter Notebook "kernel" using a new custom built environment.

2. If it's Python 2, you'll need to add another "kernel" or backend language to
   your Jupyter notebook. First we must create an "environment," which will be
   explained shortly. To do this copy the following

    ```
    conda create --yes -n cshl-sca-2017 python=3 scipy numpy jupyter scikit-learn numpy matplotlib seaborn pandas statsmodels seaborn ipykernel xlrd
    ```

3. This will create an "environment" called `cshl-sca-2017`. Environments are
   little sandboxes of software packages with different versions so they don't
   conflict with each other. You can think of adding an environment like adding
   an annex to your house where you have all the specific tools you want for a
   specific task, and all these tools are specialized and you don't want anyone
   else to touch them. For example, if you added a music recording studio to
   your house, you'd probably want headphones there, but they may be a
   completely different version of headphones than the ones you have in the
   kitchen or gym. (Sadly, my house has neither a music studio nor a gym but a
   girl can dream!)
4. **Activate the newly created environment.** This will "turn on" the
   environment you just created so we have Python 3 and all the packages we
   just added, but no more. First locate your terminal program (different for
   Mac and Windows):
   1. (Windows) After installing, search for "Anaconda" in your windows search
      bar and click on "Anaconda command prompt". In your terminal program,
      copy-paste this command and press “enter":
      ```
      activate cshl-sca-2017
      ```
   2. (Mac/Linux) Search for "Terminal" in the spotlight/program search and
      click on it. In your terminal program, copy-paste this command:
    ```
    source activate cshl-sca-2017
    ```
5. Your command prompt should now say `(cshl-sca-2017)` at the beginning of the
   line, like this:

    ```
    [obotvinnik@tscc-login2 ~]$ source activate cshl-sca-2017` discarding
    /home/obotvinnik/anaconda/bin from PATH prepending
    /home/obotvinnik/anaconda/envs/cshl-sca-2017/bin to PATH
    (cshl-sca-2017)[obotvinnik@tscc-login2 ~]$ 
    ```

6. **Add the new kernel to Jupyter Notebooks**. Now that you're in the new environment, you
   can add this Python to Jupyter as a backend language (aka "kernel") to use,
   using this command (copy the whole thing):

    ```
    python -m ipykernel install --user --display-name "Python 3 (cshl-sca-2017)"
    ```

7. **Check that it worked.** To test that the kernel added properly, start up a
   Jupyter Notebook by typing in your terminal program: `jupyter notebook`
   Which will open up a new tab in your browser. Now check that when you click
   "New," you get the option for the Python 3 kernel you just created, "Python
   3 (cshl-sca-2017)":

    ![](https://github.com/olgabot/cshl-singlecell-2017/raw/master/environment_in_jupyter_notebook_kernels.png)

<span style="background-color:cyan">STOP here if you already have Anaconda
Python installed - you are done making your environment!!</span>


Phew! That was a lot. Please email me directly (firstname.lastname at czbiohub.org) if you have any questions or issues!


Thanks,<br>
Olga