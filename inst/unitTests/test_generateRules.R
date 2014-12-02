test.generateRules <- function() {

    rsP <- new("ruleSet")
    rsP@ionlistfile <- system.file('lists/ions.csv', package = "CAMERA")
    rsP@neutraladditionfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")
    rsP@neutrallossfile <- system.file('lists/neutraladdition.csv', package = "CAMERA")

    rsP <- readLists(rsP)

    ## Empty neutraladdition, which triggers "Fehlender Wert, wo TRUE/FALSE nötig ist"
    rsP@neutraladdition <- rsP@neutraladdition[NULL,]

    rsP <- setDefaultParams(rsP)
    rsP <- generateRules(rsP)

}
