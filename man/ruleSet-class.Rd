\name{ruleSet}
\Rdversion{1.1}
\docType{class}
\alias{ruleSet-class}
\alias{class:ruleSet}

\alias{show,ruleSet-method}
\alias{generateRules}
\alias{generateRules,ruleSet-method}
\alias{setDefaultLists}
\alias{setDefaultLists,ruleSet-method}
\alias{readLists}
\alias{readLists,ruleSet-method}
\alias{setDefaultParams}
\alias{setDefaultParams,ruleSet-method}
\alias{setParams}
\alias{setParams,ruleSet,numeric,numeric,numeric,numeric,numeric,numeric,character-method}
			
\title{Class \code{ruleSet}}
\description{
  The class \code{ruleSet} is used to read lists of ions, adducts and
  neutral losses, and compile the dynamic ruleSet from those.
  This makes it possible to modify the default rules for certain
  analytical settings.
}

\section{Slots}{
  \describe{
    \item{\code{ionlistfile}:}{File of known charged ions, an example is
      found in CAMERA/lists/ions.csv .}
    \item{\code{neutrallossfile}:}{File of known neutral losses, an
      example is found in CAMERA/lists/neutralloss.csv.}    
    \item{\code{neutraladditionfile}:}{File of known adducts, an
      example is found in CAMERA/lists/lists/neutraladdition.csv .}
    \item{\code{ionlist}:}{Known charged ions.}
    \item{\code{neutralloss}:}{Known neutral losses.}
    \item{\code{neutraladdition}:}{Known adducts.}
    
    \item{\code{maxcharge}:}{.}
    \item{\code{mol	 }:}{.}
    \item{\code{nion     }:}{.}
    \item{\code{nnloss   }:}{.}
    \item{\code{nnadd    }:}{.}
    \item{\code{nh	 }:}{.}
    
    \item{\code{polarity }:}{Polarity of the ruleSet.}    
    
    \item{\code{rules}:}{data.frame of resulting mass differences, this
  is the dynamic ruleSet.}   
  }  
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  Methods implemented for \code{ruleSet}
  \describe{
    \item{setDefaultLists	}{\code{signature(object = "ruleSet")}:
     Set filenames for the lists shipped with CAMERA.}
   \item{readLists	}{\code{signature(object = "ruleSet")}: Read and
     parse the lists from the files.}
    \item{setDefaultParams}{\code{signature(object = "ruleSet")}: Set
  the default parameters for rule generation. }
    \item{setParams       }{\code{signature(object = "ruleSet")}: Set
  the parameters for rule generation.  }   
   \item{generateRules	}{\code{signature(object = "ruleSet")}: Create
  the rules in \code{ruleSet@rules} . }
  }
}
  
\author{
  Steffen Neumann and Carsten Kuhl
}

\examples{

        r <- new("ruleSet"); 
        r2 <- setDefaultLists(r) ; 
        r3 <- readLists(r2) ; 
        r4 <- setDefaultParams(r3) ; 
        r5 <- generateRules(r4)
        dim(r5@rules)

}

\keyword{classes}
