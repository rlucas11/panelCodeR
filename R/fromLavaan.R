lav2mplus <- function(lav, group.label=NULL) {

    lav_one_group <- function(lav) {
        # mplus does not like variable names with a 'dot'
        # replace them by an underscore '_'
        lav$lhs <- gsub("\\.", "_", lav$lhs)
        lav$rhs <- gsub("\\.", "_", lav$rhs)

        # remove contraints (:=, <, >, ==) here
        con.idx <- which(lav$op %in% c(":=", "<",">","=="))
        if(length(con.idx) > 0L) {
            lav <- lav[-con.idx,]
        }

        # remove exogenous variances/covariances/intercepts...
        exo.idx <- which(lav$exo == 1L & lav$op %in% c("~~", "~1"))
        if(length(exo.idx)) {
            lav <- lav[-exo.idx,]
        }

        # end of line
        lav$eol <- rep(";", length(lav$lhs))
        lav$ustart <- ifelse(is.na(lav$ustart), "", lav$ustart)
        lav$rhs2 <- ifelse(lav$free == 0L,
                           paste("@",lav$ustart,sep=""),
                           paste("*",lav$ustart,sep=""))
        lav$plabel <- gsub("\\.", "", lav$plabel)
        LABEL <- ifelse(lav$label == "", lav$plabel, lav$label)
        lav$plabel <- ifelse(LABEL == "", LABEL,
                             paste(" (", LABEL, ")",sep=""))

        # variances
        var.idx <- which(lav$op == "~~" & lav$rhs == lav$lhs)
        lav$op[var.idx] <- ""
        lav$rhs[var.idx] <- ""

        # intercepts - excluding categorical observed
        int.idx <- which(lav$op == "~1")
        lav$op[int.idx] <- ""
        lav$rhs2[int.idx] <- paste(lav$rhs2[int.idx],"]",sep="")
        lav$lhs[int.idx] <- paste("[", lav$lhs[int.idx],sep="")

        # thresholds
        th.idx <- which(lav$op == "|")
        lav$op[th.idx] <- "$"
        lav$rhs[th.idx] <- gsub("t", "", x=lav$rhs[th.idx])
        lav$rhs2[th.idx] <- paste(lav$rhs2[th.idx],"]",sep="")
        lav$lhs[th.idx] <- paste("[", lav$lhs[th.idx],sep="")

        ## Fix phantom variables
        ph.idx <- which(lav$lhs == lav$rhs & lav$op == "=~")
        lav$rhs[ph.idx] <- ""
        lav$rhs2[ph.idx] <- ""
        
        # replace binary operators
        lav$op <- ifelse(lav$op == "=~", " BY ", lav$op)
        lav$op <- ifelse(lav$op == "~", " ON ", lav$op)
        lav$op <- ifelse(lav$op == "~~", " WITH ", lav$op)

        lav2 <- paste(lav$lhs, lav$op, lav$rhs, lav$rhs2,
                      lav$plabel, lav$eol, sep="")

        body <- paste(" ", lav2, collapse="\n")

        body
    }

    body <- lav_one_group(lav)

    # constraints go to a 'MODEL CONSTRAINTS' block
    con.idx <- which(lav$op %in% c(":=", "<",">","=="))
    if(length(con.idx) > 0L) {
        ### FIXME: we need to convert the operator
        ### eg b^2 --> b**2, others??
        lav$lhs[con.idx] <- gsub("\\^","**",lav$lhs[con.idx])
        lav$rhs[con.idx] <- gsub("\\^","**",lav$rhs[con.idx])

        constraints <- "\nMODEL CONSTRAINT:\n"
        # define 'new' variables
        def.idx <- which(lav$op == ":=")
        if(length(def.idx) > 0L) {
            def <- paste(lav$lhs[def.idx], collapse= " ")
            constraints <- paste(constraints, "NEW (", def, ");")
            lav$op[def.idx] <- "="
        }
        # replace '==' by '='
        eq.idx <- which(lav$op == "==")
        if(length(eq.idx) > 0L) {
            lav$op[eq.idx] <- "="
        }
        con <- paste(gsub("\\.","",lav$lhs[con.idx]), " ",
                     lav$op[con.idx], " ",
                     gsub("\\.","",lav$rhs[con.idx]), ";", sep="")
        con2 <- paste("  ", con, collapse="\n")
        constraints <- paste(constraints, con2, sep="\n")
    } else {
        constraints <- ""
    }

    out <- paste(body, constraints, sep="")
    class(out) <- c("lavaan.character", "character")
    out
}

