# Test nohup

print(Sys.time())

for (i in 1:1000) {

  print(i)
  Sys.sleep(10)

}

# Confirm R CMD BATCH records print statements in a log.txt file: YES
# R CMD BATCH --no-save --no-restore use-nohup.R use-nohup-log.txt
