\LoadClass[a4paper,12pt]{article}

\usepackage{graphicx}
\usepackage{fancyhdr}
\usepackage{lastpage}
\usepackage{geometry}
\usepackage{adjustbox}
\usepackage{eso-pic}
\usepackage{transparent}
\usepackage{vmargin}
\usepackage{color}
\usepackage{xcolor}
\usepackage{titlesec}
\usepackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc} %Caractères spéciaux
\RequirePackage[french]{babel} %Langue du document
\usepackage[toc,page]{appendix} 
\usepackage{titletoc}
\usepackage{array}
\graphicspath{ {figures/} }
\usepackage[acronym]{glossaries}
\usepackage{float}
\usepackage{caption}
\usepackage{pdfpages}

%-------------------------------------------------------------------
%                         Mise en Page:
%-------------------------------------------------------------------
% marges
\setmarginsrb{2.5cm}{1cm}{2.5cm}{1cm}{1cm}{2cm}{1cm}{1cm}
    %1 est la marge gauche
    %2 est la marge en haut
    %3 est la marge droite
    %4 est la marge en bas
    %5 fixe la hauteur de l’entête
    %6 fixe la distance entre l’entête et le texte
    %7 fixe la hauteur du pied de page
    %8 fixe la distance entre le texte et le pied de page

%couleurs-----------------------------------------------------------
%liens
\definecolor{linkcolor}{rgb}{0.1,0.1,0.7}
\usepackage[hypertexnames=false]{hyperref}
\hypersetup{
    colorlinks,%
    citecolor=linkcolor,%
    filecolor=linkcolor,%
    linkcolor=linkcolor,%
    urlcolor=linkcolor,%
}

%titres
\definecolor{blue}{HTML}{4a86e8}
\definecolor{orange}{HTML}{ff9900}
\titleformat{\section}{\color{blue}\Large\bfseries}{\thesection}{1em}{}
\titleformat{\subsection}{\color{orange}\large\bfseries}{\thesubsection}{1em}{}



% entête
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{} % Clear header and footer
\renewcommand\headrulewidth{0.3pt}
\fancyhead[L]{\raisebox{-0.5\height}{\resizebox{!}{0.7cm}{\includegraphics{images/CRCL.png}}}} % Logo à gauche
\fancyhead[C]{\raisebox{-0.5\height}{\resizebox{!}{1.5cm}{\includegraphics{images/bioinfo.png}}}} % Logo au centre
\fancyhead[R]{\raisebox{-0.5\height}{\resizebox{!}{2cm}{\includegraphics{images/Lyon1.png}}}} % Logo à droite

% pied de page
\renewcommand\footrulewidth{0.3pt}
\fancyfoot[L]{\leftmark}

%-------------------------------------------------------------------
%                         Page de garde:
%-------------------------------------------------------------------

%--------------------------Infos------------------------------------
\newcommand{\formation}[1]{\renewcommand{\formation}{#1}}
\newcommand{\UE}[1]{\renewcommand{\UE}{#1}}
\newcommand{\titre}[1]{\renewcommand{\titre}{#1}}
\newcommand{\maitre}[1]{\renewcommand{\maitre}{#1}}
\newcommand{\tuteur}[1]{\renewcommand{\tuteur}{#1}}
\newcommand{\eleves}[1]{\renewcommand{\eleves}{#1}}
\newcommand{\annee}[1]{\renewcommand{\annee}{#1}}
%-------------------------Initialisation---------------------------
\newcommand{\pagedegarde}{
\thispagestyle{empty}
\begin{titlepage}

%wallpaper
\AddToShipoutPictureBG*{
  \AtPageLowerLeft{
    \put(\LenToUnit{\dimexpr\paperwidth-17cm},\LenToUnit{-2.5cm}){\transparent{0.3}\includegraphics[width=1\paperwidth,height=1\paperheight]{images/font.png}}
  }
}


\centering %Centraliser le contenune pas numéroter dernière page latex

%logos
\begin{figure}
\raisebox{-0.5\height}{\includegraphics[width=5cm]{images/CRCL.png}}
\hspace*{5cm}  %% pour mettre un espace (horizontal) de 5cm entre les deux images
\raisebox{-0.5\height}{\includegraphics[width=5cm]{images/Lyon1.png}}
\vspace{1cm}
\end{figure}

%titres
{\Large \formation\par}
{\normalsize \UE\par}
\vspace{1.5cm}
\begin{centering}
    \includegraphics[width=7cm]{images/bioinfo.png}
\end{centering}

\rule{\linewidth}{0.1mm}
{\huge\bfseries \titre\par}
\rule{\linewidth}{0.1mm}
\vspace{2cm}

%acteurs du projet
\begin{minipage}{0.4\textwidth}
    \begin{flushleft} \normalsize
        \emph{Maître de projet:}\par
        \maitre\par
    \end{flushleft}
    \begin{flushleft} \large
        \emph{Tuteur pédagogique:}\par
        \tuteur\par
    \end{flushleft}
\end{minipage}
\begin{minipage}{0.4\textwidth}
    \begin{flushright} \normalsize
        \emph{Étudiants:}\par
        \eleves\par
    \end{flushright}
\end{minipage}\par\vspace{2cm}

{\normalsize \annee}\par

\end{titlepage}
}

%-------------------------------------------------------------------
%                         Sommaire:
%-------------------------------------------------------------------
\newcommand{\tabledematieres}{
\newpage
\thispagestyle{fancy}
\fancyfoot[R]{}
\setcounter{page}{0}
%\tableofcontents
\startcontents
\printcontents{atoc}{0}{\section*{\contentsname}}
}

%-------------------------------------------------------------------
%                         Corps du rapport:
%-------------------------------------------------------------------
\newcommand{\corpsdetexte}{
\newpage
\thispagestyle{fancy} % Appliquer le style de page pour l'entête et le pied de page
\fancyfoot[R]{Page \thepage/\pageref{lastPage}}
\fancyfoot[L]{\textsc{\leftmark}}
}

%-------------------------------------------------------------------
%                         Bibliographie:
%-------------------------------------------------------------------
\newcommand{\bibliographie}{
\label{lastPage} %mot clé qui désigne la dernière page à numéroter
\thispagestyle{fancy}
\fancyfoot[R]{}
\bibliography{biblio}
\addcontentsline{toc}{section}{Bibliograhie}
\bibliographystyle{plain}
}

%-------------------------------------------------------------------
%                         Annexes:
%-------------------------------------------------------------------
\newcommand{\annexes}{
\newpage
\addcontentsline{toc}{section}{Annexes}
\appendix
\startcontents
\printcontents{atoc}{0}{\section*{Table des Annexes}}
\fancyfoot[L]{ANNEXES}
\newpage
\renewcommand{\thesubsection}{\Alph{subsection}}
}

%-------------------------------------------------------------------
%                         Tables:
%-------------------------------------------------------------------
\renewcommand{\listfigurename}{Liste des figures}
\renewcommand{\listtablename}{Liste des tableaux}
\newcommand{\listeTables}{
\newpage
\thispagestyle{fancy} 
\fancyfoot[R]{}
\addcontentsline{toc}{section}{Liste des tableaux}
\listoftables
}

\newcommand{\listeFigures}{
\newpage
\thispagestyle{fancy} 
\fancyfoot[R]{}
\listoffigures
\addcontentsline{toc}{section}{Liste des figures}
}

\newcommand{\listeSigles}{
\newpage
\thispagestyle{fancy} 
\fancyfoot[R]{}
\fancyfoot[L]{\MakeUppercase{Liste des abréviations, des sigles et des acronymes}}
\printglossary[type=\acronymtype]

\addcontentsline{toc}{section}{Liste des abréviations, des sigles et des acronymes}
}