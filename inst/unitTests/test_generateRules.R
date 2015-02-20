test.generateRules <- function() {

    rsP <- new("ruleSet")
    rsP@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA")
    rsP@neutraladditionfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")
    rsP@neutrallossfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")

    rsP <- readLists(rsP)

    rsP <- setDefaultParams(rsP)
    rsP <- generateRules(rsP)

    checkTrue( nrow(rsP@rules) == 89)

}


test.generateRulesNoNeutralAddition <- function() {

    rsP <- new("ruleSet")
    rsP@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA")
    rsP@neutraladditionfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")
    rsP@neutrallossfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")

    rsP <- readLists(rsP)

    ## Empty neutraladdition, which triggers "Fehlender Wert, wo TRUE/FALSE nötig ist"
    rsP@neutraladdition <- rsP@neutraladdition[NULL,]

    rsP <- setDefaultParams(rsP)
    rsP <- generateRules(rsP)

    checkTrue( nrow(rsP@rules) == 78)

}

## The following two check for
## https://github.com/sneumann/CAMERA/issues/4
##
test.generateRulesNoPosIons <- function() {
    ## Setup ruleSet
    rs <- new("ruleSet")

    rs@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA");
    rs@neutraladditionfile <- system.file('lists/neutraladdition.csv', package = "CAMERA");
    rs@neutrallossfile <- system.file('lists/neutralloss.csv', package = "CAMERA");

    rs <- readLists(rs)
    rs@ionlist<-rs@ionlist[c(2),]
    rs <- setDefaultParams(rs)

    ##default polarity is positive
    rs@polarity<-"positive"
    rs <- generateRules(rs)

    checkTrue( nrow(rs@rules) == 60)
    
}

test.generateRulesNoNegIons <- function() {
    ## Setup ruleSet
    rs <- new("ruleSet")

    rs@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA");
    rs@neutraladditionfile <- system.file('lists/neutraladdition.csv', package = "CAMERA");
    rs@neutrallossfile <- system.file('lists/neutralloss.csv', package = "CAMERA");

    rs <- readLists(rs)
    rs@ionlist<-rs@ionlist[c(1,3),]
    rs <- setDefaultParams(rs)

                                        #returns the following error in negative mode
    rs@polarity<-"negative"
    rs <- generateRules(rs)
                                        #Error in if (coeff[i, ii] > 0) { : missing value where TRUE/FALSE needed

    checkTrue( nrow(rs@rules) == 35)

}
