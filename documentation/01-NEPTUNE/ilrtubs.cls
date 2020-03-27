\NeedsTeXFormat{LaTeX2e}[1999/06/01] \ProvidesClass{ilrtubs}[2000/06/03 v0.3c IFR TU-BS]

\LoadClass[11pt,headsepline,DIV10,BCOR0.5cm]{scrbook}

%\RequirePackage{siunitx}	% SI units
\RequirePackage{amssymb}
\RequirePackage[intlimits]{amsmath}
\RequirePackage{exscale}
\RequirePackage[ansinew, latin1]{inputenc}
\RequirePackage{array}
\RequirePackage{calc}
\RequirePackage{dcolumn}
\RequirePackage{fancyhdr}
%\RequirePackage{fancyvrb}
\RequirePackage{booktabs}
\RequirePackage{dcolumn}
\RequirePackage{tabularx}
\RequirePackage{graphicx}
\RequirePackage{longtable}
%\RequirePackage{colortab}
\RequirePackage{lscape}
\RequirePackage[hang]{caption}
\RequirePackage{mathabx}
\RequirePackage{multirow}
\RequirePackage{latexsym}
\RequirePackage{textcomp}
\RequirePackage{subfigure}
%\RequirePackage[sort]{cite}
%\RequirePackage[hyphens]{url}
\RequirePackage{supertabular}
\RequirePackage{afterpage}
\RequirePackage{rotating}
\RequirePackage{flafter}
\RequirePackage{theorem}
\RequirePackage{lastpage}
\usepackage{helvet}
\usepackage{mathptmx}
\usepackage{times}
\usepackage{palatino}
\usepackage{courier}
\usepackage{pifont}

\setlength{\textheight}{42\baselineskip}
\addtolength{\textheight}{\topskip}
\newlength{\uppermargin}
\setlength{\uppermargin}{\paperheight}
\addtolength{\uppermargin}{-\textheight}
\addtolength{\uppermargin}{-\headsep}
\addtolength{\uppermargin}{-\headheight}
\setlength{\uppermargin}{.333\uppermargin} % Raender oben/unten 1:2
\addtolength{\uppermargin}{-1in}
\setlength{\topmargin}{\uppermargin}

\makeatletter

% ----- Re-definition of rules for positioning of floating objects -----
\renewcommand{\floatpagefraction}{.831}
\renewcommand{\textfraction}{.17}
\renewcommand{\topfraction}{.83}
\renewcommand{\bottomfraction}{.83}
\setcounter{topnumber}{1}                % disallow more than 1 obj. at top
\setcounter{bottomnumber}{1}             % disallow more than 1 obj. at bottom

%%% Mehr Platz nach der Zahl im Chaptertitel bei scrreprt
\renewcommand*\chapterformat {\thechapter\if@altsecnumformat.\fi\hspace{0.4cm}}

%%%----Reference anstelle von Bibliography----%%%%%%%%
\renewcommand{\bibname}{References}
%\AtBeginDocument{\renewcommand{\refname}{References}} %%%für babel, klappt aber nicht

%%%----Abbildungen nach Section und Chapter nummerieren---%%%%
%\renewcommand{\thefigure}{\thesection-\arabic{figure}}
%\@addtoreset{figure}{section}

%%%----Schriftart für Abbildungen
% caption2.sty
\renewcommand{\captionfont}{\normalfont\itshape\small}
\renewcommand{\captionlabelfont}{\normalfont\bfseries\small}

%\renewcommand{\capfont}{\itshape\small}
%\renewcommand{\caplabelfont}{\normalfont\bfseries\small}

%%%---Schraegstriche /9/ statt der eckigen Klammern [9] fuer Literatur
\renewcommand{\@cite}[2]{/{#1\if@tempswa , #2\fi}/}
\renewcommand{\@biblabel}[1]{/#1/}

%%%------Literaturverzeichnis = section statt chapter, kein Seitenwechsel--%%%%
\renewcommand*\bib@heading{%
  \section*{\bibname
    \@mkboth{\bibname}{\bibname}}}

%%%%-------Section-Schriftart für KOMA-Script ändern
%%\renewcommand*\sectfont{\sffamily}%\mdseries}

%%%-----Fig./Tab. anstatt Figure/Table-----
%\renewcommand{\figurename}{Fig.}
%\renewcommand{\tablename}{Tab.}
% --- space below caption for table captions (above the float)
\setlength{\abovecaptionskip}{3pt}
\setlength{\belowcaptionskip}{3pt}
\setlength{\textfloatsep}{10pt}
\setlength{\floatsep}{2pt}

% --- allow pagebreak within equations
%\allowdisplaybreaks[4]

%%%-------Caption-Breite von longtable---------%
\setlength{\LTcapwidth}{0.85\linewidth}
\parindent0cm       %%--keine Absatzeinrückung
\setlength{\parskip}{0.8ex}
%%--------Kürzbare Linien für Longtable--------%
   \newcommand{\LTtoprule}[1]{%
     \\
     \cmidrule[\heavyrulewidth]{#1}
   }
   \newcommand{\LTmidrule}[1]{%
     \cmidrule[\lightrulewidth]{#1}
   }
   \newcommand{\LTbottomrule}[1]{%
     \cmidrule[\heavyrulewidth]{#1}
   }
%%%-------Headings-------------%%%%
\pagestyle{fancy}

 \renewcommand{\chaptermark}[1]{%
   \markboth{\thechapter\quad#1}%
            {}%
 }

\renewcommand{\sectionmark}[1]{%   %Umstellung der Ebene die im Header angezeigt wird%
  \markright{\thesection\quad#1}%
}

\fancypagestyle{plain}{%
  \renewcommand{\headrulewidth}{0em}%
  \fancyhf{}%
  \fancyfoot[R]{\thepage}%
  \fancyfoot[L]{\small{\copyright{} Copyright TU Braunschweig, 2012}}%
}

\fancypagestyle{headings}{
  \renewcommand{\headrulewidth}{0.04em}
  \fancyhf{}
  \fancyhead[L]{\setlength{\unitlength}{1mm}
    \begin{picture}(0,0)%
      \put(0 ,-0.5){\includegraphics[height=0.9cm,clip=]{ilrlogo_englisch.eps}}%
%      \put(0 ,-0.5){\includegraphics[height=0.9cm,clip=]{ilrthumb.eps}}%
%      \put(41, 2.5){\includegraphics[height=0.25cm,clip=]{qinetiqthumb.eps}}%
%      \put(43, 5.5){\includegraphics[height=0.3cm,clip=]{aiubthumb.eps}}%
%      \put(43,-0.5){\includegraphics[height=0.18cm,clip=]{fganthumb.eps}}%
    \end{picture}
  }
% \fancyhead[LE]{\usebox{\logo}}%
%  \fancyhead[L]{\nouppercase{\footnotesize\sf{\leftmark}}}%
  \fancyhead[R]{\nouppercase{\footnotesize\sf{\rightmark}}}%
%  \fancyhead[R]{\setlength{\unitlength}{1mm}
%    \begin{picture}(0,0)%
%     \put(-18.0,-0.3){\includegraphics[width=6mm,clip=]{astrium-logo.eps}}%
%     \put(-9.5,-0.3){\includegraphics[width=1cm,clip=]{ilrlogo.eps}}%
%     \put(-18,-1){\includegraphics[width=2cm,clip=]{ifrtubs.eps}}%
%    \end{picture}
%  }
 \fancyfoot[R]{\thepage}
  \fancyfoot[L]{\small{\copyright{} Copyright TU Braunschweig, 2012}}%
}
\pagestyle{headings}

%%%%-------Änderung Paragraph-Aufbau--------
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}
  {-3.25ex\@plus -1ex \@minus -.2ex}
  {1.0ex \@plus .2ex} {\raggedsection\normalfont\size@paragraph\sectfont}}
%%%%-------Änderung Subparagraph-Aufbau--------
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\z@}
  {-2.50ex\@plus -1ex \@minus -.2ex}
  {0.1ex \@plus .2ex} {\raggedsection\scshape\size@subparagraph}}

\newcommand\refsubpar{\@startsection
  {subparagraph}{5}{\z@}%
  {-1.8ex\@plus -1ex \@minus -.2ex}%
  {0.5ex \@plus .2ex}%
  {\normalfont\normalsize\scshape}}
\newcommand\refinsubpar{\@startsection
  {subparagraph}{5}{\z@}%
  {-1.8ex\@plus -1ex \@minus -.2ex}%
  {-1.0em}%
  {\normalfont\normalsize\scshape}}

\newsavebox{\indfboxs}
\newlength{\boxstdind}\setlength{\boxstdind}{5mm}
\newlength{\boxindent}
\newlength{\indfboxwidth}
\newenvironment{indfbox}[1][\boxstdind]
  {\setlength{\indfboxwidth}{\linewidth-#1-\fboxsep-\fboxrule}%
   \setlength{\boxindent}{#1}%
   \begin{lrbox}{\indfboxs}\begin{minipage}{\indfboxwidth}}
  {\end{minipage}\end{lrbox}\noindent%
   \hspace*{\boxindent}\fbox{\usebox{\indfboxs}}}
\newcommand{\bullitm}{$\bullet$~\;}
\newlength{\bulletwidth}\settowidth{\bulletwidth}{\bullitm}
\newcommand{\bullind}{\hspace*{\bulletwidth}}

%*** definition of special itemize environment for workpackage descriptions
\newenvironment{wpditems}
{\vspace{-2mm}\begin{list} {$\bullet$}
  { \setlength{\topsep}{0pt}
    \setlength{\parskip}{0pt}
    \setlength{\partopsep}{0pt}
    \setlength{\parsep}{0pt}
    \setlength{\itemsep}{1pt}
    \setlength{\leftmargin}{4mm} }
}
{\end{list}\vspace{-5mm}}
\newenvironment{wpditemx}
{\vspace{0mm}\begin{list} {--}
  { \setlength{\topsep}{0pt}
    \setlength{\parskip}{0pt}
    \setlength{\partopsep}{0pt}
    \setlength{\parsep}{0pt}
    \setlength{\itemsep}{1pt}
    \setlength{\leftmargin}{4mm} }
}
{\end{list}\vspace{0mm}}
\newenvironment{cvitems}
{\vspace{-3mm}\begin{list} {$\bullet$}
  { \setlength{\topsep}{0pt}
    \setlength{\parskip}{0pt}
    \setlength{\partopsep}{0pt}
    \setlength{\parsep}{0pt}
    \setlength{\itemsep}{1pt}
    \setlength{\leftmargin}{4mm} }
}
{\end{list}}
\newcommand{\cvskip}{\\[-4.2mm]}
\newcommand{\wpdtail}[2]{
\begin{tabular}{|p{23.4mm}||p{64.3mm}||p{40mm}|}
  \textbf{\rule[-1mm]{0mm}{5mm}WP Manager} & #1 & 
     \textbf{Effort}: #2\\
  \hline\hline
\end{tabular}
\vspace{-0.5mm}
}
\newcommand{\wpditem}[3]{
\begin{tabular}{|p{138.2mm}|}
  \rule{0mm}{5mm}\textbf{#1:}\\
  #2\\[1mm]
  #3
\end{tabular}
\vspace{-0.5mm}
}
\newcommand{\wpdobjc}[1]{
  \wpditem{Objectives}{#1}{}
}
\newcommand{\wpdinpt}[1]{
  \wpditem{Required Input}{#1}{}
}
\newcommand{\wpdintf}[1]{
  \wpditem{Restraining Interfaces}{#1}{}
}
\newcommand{\wpdtask}[1]{
  \wpditem{Tasks}{#1}{}
}
\newcommand{\wpdoutp}[1]{
  \wpditem{Required Output}{#1}{\hline}
}
\newcommand{\cwp}{tbd.}

%*** definition of special itemize environment for curriculae vitae
\newcommand{\cvitem}[4]{
\renewcommand{\arraystretch}{1.2}
\begin{tabular}{p{24mm}p{110mm}}
  #1  #2
  #3\\#4
\end{tabular}
}
\newcommand{\cvname}[1]{
  \cvitem{}{\rule[-5mm]{0mm}{5mm}\textbf{Name:} & \textbf{#1}}{}{\hline}
}
\newcommand{\cvborn}[1]{
  \cvitem{}{\rule{0mm}{7mm}\textbf{Date of Birth:} & #1}{}{}
}
\newcommand{\cveduc}[1]{
  \cvitem{}{\textbf{Qualification:} & \\}{#1}{}
}
\newcommand{\cvlang}[1]{
  \cvitem{}{\textbf{Languages:} & #1}{}{}
}
\newcommand{\cvassg}[1]{
  \cvitem{}{\multicolumn{2}{l}{\textbf{Project Assignment:}} \\}{& #1}{}
}
\newcommand{\cvprof}[1]{
  \cvitem{}{\multicolumn{2}{l}{\textbf{Professional History:}} \\}{#1}{}
}
\newcommand{\cvproj}[1]{
  \cvitem{}{\multicolumn{2}{l}{\textbf{Project History:}} \\}{#1}{}
}

%%%--- Anpassung der Einstellungen fuer schoenere Tabellen
\newcommand{\armultirow}[3]{%
  \multicolumn{#1}{#2}{%
    \begin{picture}(0,0)%
      \put(0,0){%
        \begin{tabular}[t]{@{}#2@{}}%
          #3%
        \end{tabular}%
      }%
    \end{picture}%
  }%
}%
%
% array.sty
\newcolumntype{f}[1]{>{\raggedright\centering\hspace{0pt}}p{#1}}
\newcolumntype{g}{>{$}r<{$}}
\newcolumntype{n}{l}
\newcolumntype{N}{>{\scriptsize}l}
\newcolumntype{v}[1]{>{\raggedright\hspace{0pt}}p{#1}}
\newcolumntype{u}{>{\ttfamily}l}
\newcolumntype{s}{>{\sf}l}
\newcolumntype{V}[1]{>{\scriptsize\raggedright\hspace{0pt}}p{#1}}

%dcolum.sty
%\newcolumntype{d}{D{ d}{ d}{-1}}
\newcolumntype{d}[1]{>{\DC@{.}{.}{#1}}c<{\DC@end}}

% array.sty, tabularx.sty
\newcolumntype{x}{>{\scriptsize\raggedright\hspace{0pt}}X}
%
% array.sty, rotating.sty
\newcolumntype{R}[1]{%
  >{\begin{turn}{90}\begin{minipage}{#1}\raggedright\hspace{0pt}}c%
  <{\end{minipage}\end{turn}}%
}

%%%%---"Entspricht"-Zeichen----
\newcommand{\entspricht}{\mathrel{\widehat{=}}}

%clear page for new chapter on odd page
\newcommand{\chapnew}{\thispagestyle{plain}\cleardoublepage}

%short cuts for major software modules
%\newcommand{\sub}[1]{\textsc{#1}}
%\newcommand{\file}[1]{\texttt{\textit{#1}}}

\makeatother
\setcounter{secnumdepth}{3}  %% Numerierungsebenen: 1.1.4. Wert erhoehen, wenn tiefer...(article)
\setcounter{tocdepth}{3}     %% dto. fuer Inhaltsverzeichnis, wieviele Ebenen aufgenommen werden
\pagenumbering{roman}        %% roemische Zahlen i, ii, iii fuers Inhaltsverzeichnis, {Roman} fuer grosse Zahlen
\newpage
\frenchspacing


% ----- short cuts for certain file types -----
\newcommand{\program}[1]{\emph{#1}}
\newcommand{\sub}[1]{\textsc{#1}}
\newcommand{\file}[1]{\texttt{\textit{#1}}}
\newcommand{\auth}[1]{\textsl{#1}}
\newcommand{\var}[1]{\emph{#1}}               % variables
\newcommand{\outp}[1]{\texttt{#1}}            % program input/output
\newcommand{\tit}[1]{\emph{#1}}               % document titles
\newcommand{\eqv}{\mathrel{\widehat{=}}}
\newcommand{\eff}{\text{eff}}
\newcommand{\erf}{\text{erf}}

% ----- short cuts for major software modules -----
\newcommand{\master}{MASTER\/}
\newcommand{\masterd}{\master{}-2001\/}
\newcommand{\mastere}{\master{}-2005\/}
\newcommand{\masterf}{\master{}-2009\/}
\newcommand{\mastero}{\masterd{}}
\newcommand{\mastern}{\masterf{}}
%\newcommand{\mastera}{\master{}$^{2009}$\/}
\newcommand{\analyst}{\program{ANALYST}\/}
\newcommand{\standard}{\program{STANDARD}\/}
\newcommand{\poem}{\program{POEM}\/}
\newcommand{\poemd}{\poem{}-2005\/}
\newcommand{\proof}{\program{PROOF}\/}
\newcommand{\proofd}{\proof{}-2001\/}
\newcommand{\proofe}{\proof{}-2005\/}
\newcommand{\prooff}{\proof{}-2009\/}
\newcommand{\proofn}{\prooff{}}
%\newcommand{\proofa}{\proof{}$^{2009}$\/}
\newcommand{\focus}{\program{FOCUS1}\/}
\newcommand{\metmod}{\program{METMOD}\/}
\newcommand{\disco}{\program{DISCO}\/}
\newcommand{\probdens}{\program{PROBDENS}\/}
\newcommand{\qdelta}{\program{DELTA}\/}
\newcommand{\deltaqq}{\qdelta\/}
\newcommand{\gnuplot}{\textbf{gnuplot}\/}
\newcommand{\fluxbrowser}{\program{Flux Browser}\/}
\newcommand{\marwin}{\program{MARWIN}\/}

% ---- short cut for MASTER reference epoch
\newcommand{\tref}{May\,1,\:2009}

% ---- colored cell commands
\newcommand{\cocell}[4]{%
%\LCC
%\lightgray & \lightgray & \lightgray & \lightgray \\
\rule[-3mm]{0mm}{8mm}\textbf{#1} & \textbf{#2} & \textbf{#3} & \textbf{#4}\\ \RP \hline
%\ECC
}

% ---- generic imagewidth
\newlength{\imagewidth}
\setlength{\imagewidth}{0.7\textwidth}
\newlength{\smallimagewidth}
\setlength{\smallimagewidth}{0.6\textwidth}


% ---- additional math operators
\DeclareMathOperator{\artanh}{artanh}
\DeclareMathOperator{\arcosh}{arcosh}
\newcommand{\tm}{{\texttrademark}}
\newcommand{\rr}{{\textregistered}}


% ---- generic verbatim text style
%\newenvironment{verbfoot}{\footnotesize\verbatim}
%                         {\endverbatim}
%\newenvironment{verbscript}{\scriptsize\verbatim}
%                           {\endverbatim\vspace*{1ex}}
