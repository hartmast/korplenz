
# Zusatzpakete installieren und laden
sapply(c("dplyr", "devtools"), 
       function(x) if(!is.element(x, installed.packages())) 
         install.packages(x, dependencies = T))

# collostructions
if(!is.element("collostructions", installed.packages())) {
  if(.Platform$OS.type!="Windows") {
    install.packages("http://userpage.fu-berlin.de/~flach/wp-content/uploads/collostructions_0.0.10.tar.gz", repos = NULL)
  } else {
    install.packages("http://userpage.fu-berlin.de/~flach/wp-content/uploads/collostructions_0.0.10.zip", repos = NULL)
  }
}

# concordances
if(!is.element("concordances", installed.packages())) {
  devtools::install_github("hartmast/concordances")
}

library(dplyr)
library(concordances)
library(collostructions)

# Optionen
Sys.setlocale("LC_ALL", "de_DE") # funktioniert nur bei Mac/Linux und eruebrigt sich meist
options(stringsAsFactors = F)

# Daten einlesen und aufbereiten:

# Konkordanz
allein <- getCOSMAS("selltenallein_WoN.txt.TXT", years = T)

# Wortformenliste: "DeReKo-2014-II-MainArchive-STT.100000.freq"
# verfuegbar hier: http://www1.ids-mannheim.de/kl/projekte/methoden/derewo.html
wf <- read.csv(file.choose(), sep="\t", head=F,
               encoding="Latin-1", quote="")
colnames(wf) <- c("Token", "Lemma", "POS", "Freq")

# Interpunktion entfernen
allein$Key2 <- gsub("[[:punct:]]", "", allein$Key2)

# Tabelle mit Frequenzen der Nomina im N-Slot von [ein N kommt selten allein]
tbl <- allein$Key2 %>% table %>% sort(decreasing=T) %>% as.data.frame
colnames(tbl) <- c("word", "cxn.freq")

# Gesamtfrequenz der einzelnen Woerter im Korpus hinzufügen

# neue Spalte für die Gesamtfrequenz
tbl$word.freq <- NA

# neue Spalte auffuellen
for(i in 1:nrow(tbl)) {
  if(tbl$word[i] %in% filter(wf, POS %in% c("NN", "NE"))$Token) {
    # Lemma finden
    la <- wf[which(wf$Token==tbl$word[i] & wf$POS %in% c("NN", "NE")),]$Lemma[1]
    
    # Frequenz des Lemmas ermitteln
    tbl$word.freq[i] <- sum(wf[which(wf$Lemma==la & wf$POS %in% c("NN", "NE")),]$Freq)
  }
}

# NAs entfernen
tbl <- tbl[!is.na(tbl$word.freq),]

# "corpsize" ermitteln (s. ?collex): hier, Gesamtfrequenz aller Substantive
corpsize <- sum(filter(wf, POS %in% c("NN", "NE"))$Freq)

# collex
collex(tbl, corpsize = corpsize) %>% write.table("collex_results.csv",
                                                 sep="\t", row.names=F,
                                                 quote=F,
                                                 fileEncoding="UTF-8")




