## Script to analyze the risk/warning levels produced by the forecasting, and judge if a trigger should be created to distribute the information.
## By Jafet Andersson, Feb, 2019
7
# The logic of the script
# 1. Analyze the "...forecast_mapWarningLevel.txt" file
# 2. Apply a few criteria for wether a trigger should be made or not
# 3. Create the trigger & and the messages to be sent by SMS and EMAIL.


# To do
# - add dynamic info for the model (now hardcoded to NH_2.23, but in future take it from some log information. BRIEF)
# - Add link to the job (Emmanuel?), read from some input par?
# - make nicers templates (images for email, insert all text from the criterion dynamically to not have so many templates?)
# - different settings (e.g. different criteria, different regions, different send lists)

download_file <- function(local_file_name, url) {

  # delete local file if it existst
  if(file.exists(local_file_name)){
     file.remove(local_file_name)
  }

  # Copy apikey if not exist
  if(!file.exists(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "apikey", sep="/"))){
    rciop.log ("ERROR", paste0("Please create the file 'apikey' in _CIOP_APPLICATION_PATH root, usually at ", Sys.getenv("_CIOP_APPLICATION_PATH"), "/apikey"))
    quit(status=1)
  }

  f_con <- file(paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "apikey", sep="/"), "r")
  first_line <- readLines(f_con,n=1)
  close(f_con)

  
  apikey <- trimws(first_line)

  # system command to download remote file to local path
  sysCmd=paste("curl", "-u", apikey, "-o", local_file_name,url,sep=" ")
  system(sysCmd,intern=T)  

  return (local_file_name)
}

TriggerDistribution <- function(path, idate) {
    
    # Define output list
    output_files <- list() 

    # --------------------------------
    # 1 - General
    path.wl<-paste(path, "/", sep="")  # define the path to the location of  "...forecast_mapWarningLevel.txt"
    name.wl<-dir(path.wl,pattern= "forecast_mapWarningLevel.txt")

    name.wl_out <- gsub(substr(name.wl, 1, 8), as.character(idate), name.wl, fixed=T)
    name.wl_out <- paste0(substr(name.wl, 1, 4), name.wl_out)

    # Read this from the hydro-smhi fileshare store
    #path.templates <- paste(Sys.getenv("_CIOP_APPLICATION_PATH"), "util/R/triggercode/templates/", sep="/")
    path.templates <- "./"
    
    download_file("./sms_nigerhype.txt", "'https://store.terradue.com/hydro-smhi/fanfar/distribution-templates/sms_nigerhype.txt'")
    download_file("./email_nigerhype.txt", "'https://store.terradue.com/hydro-smhi/fanfar/distribution-templates/email_nigerhype.txt'")

    #sysCmd=paste("mv", "./sms_nigerhype.txt", path.templates, sep=" ")
    #sysCmd=paste("mv", "./email_nigerhype.txt", path.templates, sep=" ")



    #path.templates<-"../util/R/triggercode/templates/"  # path to where the message templates are stored. On our Store
    name.temp.sms<-"sms_nigerhype.txt"  # name of the SMS template, make dynamic, reading an input parameter instead for which model to use
    name.temp.email<-"email_nigerhype.txt"  # dito for email, possibly change to some format that can read images?
    path.out <- paste0(path.wl, name.wl_out)  # path to where the outputs shall be written. Normally the same as path.wl

    
    # ---------------------------------
    # 2 - Process
    # get data
    wl<-read.table(paste0(path.wl,name.wl),header=T,sep=",",skip=1)
    temp.sms<-readLines(paste0(path.templates,name.temp.sms))
    temp.email<-readLines(paste0(path.templates,name.temp.email))
    
    # determine how many subbasins are at or above a certain warning level
    wl1<-length(which(wl$WarningLevel>=1))
    wl2<-length(which(wl$WarningLevel>=2))
    wl3<-length(which(wl$WarningLevel>=3))
    

    # trigger functions
    # If >10% of subbasins are at or above warning level 2
    t2fun<-function(mywl2) { 
        myval<-mywl2/nrow(wl)
        if(myval>0.1) TRUE else FALSE
        }
    # If any subbasins are are at or above Warning level 3
    t3fun<-function(mywl3) { 
        if(mywl3>0) TRUE else FALSE
        }
    
    
    # Run the trigger functions & prepare template
    if (t2fun(wl2) | t3fun(wl3)) {  # any trigger yes/no
        # initiate messages
        sms<-temp.sms[1]
        email<-temp.email[1:2]
        email[2]<-sub("YYYYMMDD",paste(unlist(strsplit(name.wl,split="_"))[2:3],collapse=" "),email[2])  # insert idate instead
        
        if (t2fun(wl2)) { # trigger t2fun
            sms[length(sms)+1]<-sub("XX",wl2,temp.sms[2])  # todo: don't hardcode the positions of the different triggers
            email[length(email)+1]<-sub("XX",wl2,temp.email[3]) # dito
            }
        if (t3fun(wl3)) { # trigger t3fun, todo: don't hardcode the positions of the different triggers
            sms[length(sms)+1]<-sub("YY",wl3,temp.sms[3])
            email[length(email)+1]<-sub("YY",wl3,temp.email[4])
            }
        
        # Add link here # fixme
        sms<-c(sms,temp.sms[4])
        email<-c(email,temp.email[5:6])
        
        # write trigger file and message files, change to idate in the file names
        write("Triggers activated",file=paste0(path.out,"_send_messages.txt"))  # fixme: make dynamic to depend on the model/job instead
        writeLines(text=sms,con=paste0(path.out,"_sms_message.txt"))
        writeLines(text=email,con=paste0(path.out,"_email_message.txt"))
        
        output_files <- c(output_files, paste0(path.out,"_send_messages.txt"))
        output_files <- c(output_files, paste0(path.out,"_sms_message.txt"))
        output_files <- c(output_files, paste0(path.out, "_email_message.txt"))

        #output_files <- c(output_files, paste0(path.out,paste0(paste(unlist(strsplit(name.wl,split="_"))[2:3],collapse="_"),"_email_message.txt")))

    } else {
        # write non-trigger file (do it since then one can always check if the process worked or not)  
        write("No triggers activated",file=paste0(path.out,"_donotsend_messages.txt"))  # fixme: make dynamic to depend on the model/job instead
        output_files <- c(output_files, paste0(path.out,"_donotsend_messages.txt"))
    }
        return (output_files)
}
  
  
