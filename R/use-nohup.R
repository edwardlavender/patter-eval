# Test nohup

for (i in 1:1000) {

  print(i)

}

# Confirm R CMD BATCH records print statements in a log.txt file: YES
# R CMD BATCH --no-save --no-restore use-nohup.R use-nohup-log.txt
