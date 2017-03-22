% -------------------------------------------------------------------
%  @LaTeX-class-file{
%     filename        = "Dissertate.cls",
%     version         = "2.0",
%     date            = "25 March 2014",
%     codetable       = "ISO/ASCII",
%     keywords        = "LaTeX, Dissertate",
%     supported       = "Send email to suchow@post.harvard.edu.",
%     docstring       = "Class for a dissertation."
% --------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Dissertate}[2014/03/25 v2.0 Dissertate Class]
\LoadClass[12pt, oneside, letterpaper]{book}

% Text layout.
\RequirePackage[width=5.75in, letterpaper]{geometry}
\usepackage{ragged2e}
\RaggedRight
\RequirePackage{graphicx}
\usepackage{fixltx2e}
\parindent 12pt
\RequirePackage{lettrine}
\RequirePackage{setspace}
\RequirePackage{verbatim}

\usepackage{upgreek} % SJF, 3/16/17
\usepackage{tabularx} %Edlyn added this, Jan 28
\newcolumntype{Y}{>{\centering\arraybackslash}X}%Edlyn added this, Jan 28
\newcolumntype{b}{X}%Edlyn added this, Jan 28
\newcolumntype{s}{>{\hsize=.45\hsize}X}%Edlyn added this, Jan 28

\usepackage{enumitem} %Edlyn added this, March 23
\usepackage{indentfirst} % SJF

% Fonts.
\RequirePackage{color}
\usepackage{hyperref}
\RequirePackage{url}
\RequirePackage{amssymb}
\usepackage{amsmath} %Edlyn Jan 12
\RequirePackage{mathspec}
\setmathsfont(Digits,Latin,Greek)[Numbers={Proportional}]{TeX Gyre Termes}
\setmathrm{TeX Gyre Termes}
\widowpenalty=300
\clubpenalty=300
\setromanfont[Numbers=OldStyle, Ligatures={Common, TeX}, Scale=1.0]{TeX Gyre Termes}
\newfontfamily{\smallcaps}[RawFeature={+c2sc,+scmp}]{TeX Gyre Termes}
\setsansfont[Scale=MatchLowercase, BoldFont={Lato Bold}]{Lato Regular}
\setmonofont[Scale=MatchLowercase]{Source Code Pro}
\RequirePackage[labelfont={bf,sf,footnotesize,singlespacing},
                textfont={sf,footnotesize,singlespacing},
                justification={justified,RaggedRight},
                singlelinecheck=false,
                margin=0pt,
                figurewithin=chapter,
                tablewithin=chapter]{caption}
\renewcommand{\thefootnote}{\fnsymbol{footnote}}

% Headings and headers.
\RequirePackage{fancyhdr}
\RequirePackage[tiny, md, sc]{titlesec}
\setlength{\headheight}{15pt}
\pagestyle{plain}
\RequirePackage{titling}

% Front matter.
\setcounter{tocdepth}{1}
\usepackage[titles]{tocloft}
\usepackage[titletoc]{appendix}
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}
\renewcommand{\cftchapfont}{\normalsize \scshape}
\renewcommand\listfigurename{List of figures}
\renewcommand\listtablename{List of tables}

% Endmatter
\renewcommand{\setthesection}{\arabic{chapter}.A\arabic{section}}

% References.
\renewcommand\bibname{References}
\RequirePackage[super,comma,numbers,sort&compress]{natbib} %"sort&compress" Added by Edlyn, Nov. 29 -- cv to sort and compress multiple citation numbers 
% \RequirePackage[square,comma,numbers,sort&compress]{natbib} % SJF
\renewcommand{\bibnumfmt}[1]{[#1]}
\RequirePackage[palatino]{quotchap}
\renewcommand*{\chapterheadstartvskip}{\vspace*{-0.5\baselineskip}}
\renewcommand*{\chapterheadendvskip}{\vspace{1.3\baselineskip}}

% An environment for paragraph-style section.
\providecommand\newthought[1]{%
   \addvspace{1.0\baselineskip plus 0.5ex minus 0.2ex}%
   \noindent\textsc{#1}%
}

% Align reference numbers so that they do not cause an indent.
\newlength\mybibindent
\setlength\mybibindent{0pt}
\renewenvironment{thebibliography}[1]
    {\chapter*{\bibname}%
     \@mkboth{\MakeUppercase\bibname}{\MakeUppercase\bibname}%
     \list{\@biblabel{\@arabic\c@enumiv}}
          {\settowidth\labelwidth{\@biblabel{999}}
           \leftmargin\labelwidth
            \advance\leftmargin\dimexpr\labelsep+\mybibindent\relax\itemindent-\mybibindent
           \@openbib@code
           \usecounter{enumiv}
           \let\p@enumiv\@empty
           \renewcommand\theenumiv{\@arabic\c@enumiv}}
     \sloppy
     \clubpenalty4000
     \@clubpenalty \clubpenalty
     \widowpenalty4000%
     \sfcode`\.\@m}
    {\def\@noitemerr
      {\@latex@warning{Empty `thebibliography' environment}}
     \endlist}

% Some definitions.
\def\degreeyear#1{\gdef\@degreeyear{#1}}
\def\degreemonth#1{\gdef\@degreemonth{#1}}
\def\degree#1{\gdef\@degree{#1}}
\def\advisor#1{\gdef\@advisor{#1}}
\def\department#1{\gdef\@department{#1}}
\def\field#1{\gdef\@field{#1}}
\def\university#1{\gdef\@university{#1}}
\def\universitycity#1{\gdef\@universitycity{#1}}
\def\universitystate#1{\gdef\@universitystate{#1}}
\def\programname#1{\gdef\@programname{#1}}
\def\pdOneName#1{\gdef\@pdOneName{#1}}
\def\pdOneSchool#1{\gdef\@pdOneSchool{#1}}
\def\pdOneYear#1{\gdef\@pdOneYear{#1}}
\def\pdTwoName#1{\gdef\@pdTwoName{#1}}
\def\pdTwoSchool#1{\gdef\@pdTwoSchool{#1}}
\def\pdTwoYear#1{\gdef\@pdTwoYear{#1}} % sjf
\def\pdThreeName#1{\gdef\@pdThreeName{#1}} % sjf
\def\pdThreeSchool#1{\gdef\@pdThreeSchool{#1}} % sjf
\def\pdThreeYear#1{\gdef\@pdThreeYear{#1}} % sjf
\def\pdFourName#1{\gdef\@pdFourName{#1}} % sjf
\def\pdFourSchool#1{\gdef\@pdFourSchool{#1}} % sjf
\def\pdFourYear#1{\gdef\@pdFourYear{#1}} % sjf
% School name and location
\university{Harvard University}
\universitycity{Cambridge}
\universitystate{Massachusetts}

% School color found from university's graphic identity site:
% http://isites.harvard.edu/icb/icb.do?keyword=k75408&pageid=icb.page392732
\definecolor{SchoolColor}{rgb}{0.6471, 0.1098, 0.1882} % Crimson
\definecolor{chaptergrey}{rgb}{0.6471, 0.1098, 0.1882} % for chapter numbers

\hypersetup{
    colorlinks,
    citecolor=SchoolColor,
    filecolor=black,
    linkcolor=black,
    urlcolor=SchoolColor,
}

% Formatting guidelines found in:
% http://www.gsas.harvard.edu/publications/form_of_the_phd_dissertation.php
\renewcommand{\frontmatter}{
	\input{frontmatter/personalize}
	\maketitle
	\copyrightpage
	\abstractpage
	\tableofcontents
	\listoffigures % optional
	\dedicationpage
	\acknowledgments
}

\renewcommand{\maketitle}{
	\thispagestyle{empty}
	\vspace*{\fill}
	\vspace{100pt}
	\begin{center}
	\Huge \textcolor{SchoolColor}{\thetitle} \normalsize \\
	\vspace{100pt}
	\textsc{a dissertation presented \\ by\\
	\theauthor\\ to\\ The Department of \@department\\
	\vspace{12pt}
	in partial fulfillment of the requirements\\
	for the degree of\\ 
	\@degree\\
	in the subject of\\ 
	\vspace{12pt}
	\@field\\
	\vspace{12pt}
	\@university\\ \@universitycity, \@universitystate\\
	\@degreemonth\ \@degreeyear}
	\end{center} \vspace*{\fill}
}

\newcommand{\copyrightpage}{
	\newpage
	\thispagestyle{empty}
	\vspace*{\fill}
	\scshape \noindent \small \copyright \small \@degreeyear \hspace{3pt}-- \theauthor \\
	\noindent all rights reserved.
	\vspace*{\fill}
	\newpage
	\rm
}

\newcommand{\abstractpage}{
	\newpage
	\pagenumbering{roman}
	\setcounter{page}{3}
	\addcontentsline{toc}{chapter}{Abstract} %Edlyn added Nov 18, 2015
	\pagestyle{fancy}
	\lhead{Thesis advisor: Professor \@advisor} \rhead{Author: \@author} %Edlyn added "Author:" Nov 17, 2015
	\renewcommand{\headrulewidth}{0.0pt}
	\vspace*{35pt}
	\begin{center}
	\Large \textcolor{SchoolColor}{\@title} \normalsize \\
	\vspace*{20pt}
	\scshape Abstract \\ \rm
	\end{center}
	\input{frontmatter/abstract}
	\vspace*{\fill}
	\newpage \lhead{} \rhead{}
	\cfoot{\thepage}
}

\newcommand{\dedicationpage}{
	\newpage \thispagestyle{fancy} \vspace*{\fill}
	\addcontentsline{toc}{chapter}{Dedication} %Edlyn added Nov 18, 2015
	\scshape \noindent \input{frontmatter/dedication}
	\vspace*{\fill} \newpage \rm
}

\newcommand{\acknowledgments}{
	\chapter*{Acknowledgments}
	\noindent
	\addcontentsline{toc}{chapter}{Acknowledgments} %Edlyn added Nov 18, 2015
	\input{frontmatter/acknowledgments}
	\vspace*{\fill} \newpage
	\setcounter{page}{1}
	\pagenumbering{arabic}
}


%These are new commands added by Edlyn for mathematics formatting 
\newcommand{\diff}{\mathrm{d}} %upright d for derivatives
\newcommand\fulld[2]{\frac{\mathrm{d} #1}{\mathrm{d} #2}} %full derivative
\newcommand\pard[2]{\frac{\partial #1}{\partial #2}} %partial derivative
\newcommand\pardd[2]{\frac{\partial^2 #1}{\partial #2^2}} %partial derivative
\newcommand{\vgrad}{\pmb{\nabla}} %grad
\newcommand{\vdiv}{\pmb{\nabla}\cdot} %divergence
\newcommand{\vcurl}{\pmb{\nabla}\times} %curl
\newcommand\tensor[1]{\overline{\overline{#1}}} %double bar tensor notation

% SJF
\newcommand{\unit}[1]{\textsf{#1}} % font formatting for units
\newcommand{\dna}[1]{\texttt{#1}} % font formatting for dna sequence