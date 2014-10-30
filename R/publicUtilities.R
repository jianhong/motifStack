readPCM <- function(path=".", pattern=NULL){
    pcms <- dir(path,pattern)
    pcml <- lapply(pcms, function(.ele){
        data <- read.table(file.path(path, basename(.ele)))
        classes <- sapply(data, class)
        data <- data[, classes %in% c("integer", "numeric")]
        rownames(data) <- c("A", "C", "G", "T")
        data
    })
    names(pcml) <- gsub("\\.(pcm|txt)$", "", basename(pcms), ignore.case=TRUE)
    pcm <- mapply(function(.d, .n){
        new("pcm", mat=as.matrix(.d), name=gsub("\\.", "_", make.names(.n))) 
    }, pcml, names(pcml))
    names(pcm) <- gsub("\\.", "_", make.names(names(pcm)))
    pcm
}

readPWM <- function(path=".", pattern=NULL, to=c("pfm", "pcm")){
    to <- match.arg(to)
    pwms <- dir(path,pattern)
    pwml <- lapply(pwms, function(.ele){
        data <- readLines(file.path(path, basename(.ele)))
        data <- data[grepl("^[ACGT]:",data)]
        data <- do.call(rbind, strsplit(data, "\\t"))
        rownames(data) <- gsub(":", "", data[,1])
        data <- data[c("A","C","G","T"), -1]
        data <- as.data.frame(data)
        for(i in 1:ncol(data)) data[,i] <- as.numeric(as.character(data[,i]))
        as.matrix(data)
    })
    names(pwml) <- gsub("[\\._](pwm|txt|pwm\\.txt)$", "", basename(pwms), ignore.case=TRUE)
    pwm <- mapply(function(.d, .n){
        if(to=="pfm") new("pfm", mat=as.matrix(.d), name=gsub("\\.", "_", make.names(.n)))
        else{
            .d <- round(.d * 1000)
            new("pcm", mat=as.matrix(.d), name=gsub("\\.", "_", make.names(.n)))
        }
    }, pwml, names(pwml))
    names(pwm) <- gsub("\\.", "_", make.names(names(pwm)))
    pwm
}

colorset<-function(alphabet="DNA", colorScheme='auto'){
    if(!alphabet %in% c("DNA","RNA","AA")) stop("alphabet must be one of 'DNA', 'RNA' or 'AA'")
    if(alphabet=='PROTEIN' & !(colorScheme %in% c('auto', 'charge', 'chemistry', 'classic', 'hydrophobicity')))
    stop("color scheme must be one of 'auto', 'charge', 'chemistry', 'classic' or 'hydrophobicity' for protein")
    if(alphabet %in% c('DNA','RNA') & !(colorScheme %in% c('auto', 'basepairing')))
    stop("color scheme must be one of 'auto' or 'basepairing'")
    taylor<-c(  'A'='#CCFF00',
                'C'='#FFFF00',
                'D'='#FF0000',
                'E'='#FF0066',
                'F'='#00FF66',
                'G'='#FF9900',
                'H'='#0066FF',
                'I'='#66FF00',
                'K'='#6600FF',
                'L'='#33FF00',
                'M'='#00FF00',
                'N'='#CC00FF',
                'P'='#FFCC00',
                'Q'='#FF00CC',
                'R'='#0000FF',
                'S'='#FF3300',
                'T'='#FF6600',
                'V'='#99FF00',
                'W'='#00CCFF',
                'Y'='#00FFCC')
    charge<-c(   'A'='#CCCCCC',
                 'C'='#CCCCCC',
                 'D'='#FFB32C',
                 'E'='#FFB32C',
                 'F'='#CCCCCC',
                 'G'='#CCCCCC',
                 'H'='#2000C7',
                 'I'='#CCCCCC',
                 'K'='#2000C7',
                 'L'='#CCCCCC',
                 'M'='#CCCCCC',
                 'N'='#CCCCCC',
                 'P'='#CCCCCC',
                 'Q'='#CCCCCC',
                 'R'='#2000C7',
                 'S'='#CCCCCC',
                 'T'='#CCCCCC',
                 'V'='#CCCCCC',
                 'W'='#CCCCCC',
                 'Y'='#CCCCCC')
    chemistry<-c(   'A'='#000000',
                    'C'='#00811B',
                    'D'='#D00001',
                    'E'='#D00001',
                    'F'='#000000',
                    'G'='#00811B',
                    'H'='#2000C7',
                    'I'='#000000',
                    'K'='#2000C7',
                    'L'='#000000',
                    'M'='#000000',
                    'N'='#800080',
                    'P'='#000000',
                    'Q'='#800080',
                    'R'='#2000C7',
                    'S'='#00811B',
                    'T'='#00811B',
                    'V'='#000000',
                    'W'='#000000',
                    'Y'='#00811B')
    hydrophobicity<-c(   'A'='#00811B',
                         'C'='#2000C7',
                         'D'='#000000',
                         'E'='#000000',
                         'F'='#2000C7',
                         'G'='#00811B',
                         'H'='#00811B',
                         'I'='#2000C7',
                         'K'='#000000',
                         'L'='#2000C7',
                         'M'='#2000C7',
                         'N'='#000000',
                         'P'='#00811B',
                         'Q'='#000000',
                         'R'='#000000',
                         'S'='#00811B',
                         'T'='#00811B',
                         'V'='#2000C7',
                         'W'='#2000C7',
                         'Y'='#2000C7')
    base_pairing<-c('A'="#ff8c00",'C'="#2000C7",'G'="#2000C7",'TU'="#ff8c00")
    nucleotide<-c('A'="#00811B",'C'="#2000C7",'G'="#FFB32C",'TU'="#D00001")
    if(alphabet=='DNA'){
        color<-switch(colorScheme, auto=nucleotide, basepairing=base_pairing)
        names(color)<-c('A','C','G','T')
        color
    }else{
        if(alphabet=='RNA'){
            color<-switch(colorScheme, auto=nucleotide, basepairing=base_pairing)
            names(color)<-c('A','C','G','U')
            color
        }else{
            switch(colorScheme, 
                   auto=taylor, 
                   charge=charge, 
                   chemistry=chemistry, 
                   classic=taylor, 
                   hydrophobicity=hydrophobicity)
        }
    }
}


highlightCol <- function (col, alpha=0.5){
    n <- names(col)
    col <- grDevices::col2rgb(col,alpha=T)
    col <- apply(col, 2, function(.ele, alpha){rgb(.ele[1], .ele[2], .ele[3], alpha=ceiling(alpha*.ele[4]), maxColorValue=255)}, alpha)
    col <- unlist(col)
    names(col) <- n
    col
}