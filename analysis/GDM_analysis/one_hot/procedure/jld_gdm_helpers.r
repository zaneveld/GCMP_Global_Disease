## GDM helper functions
    # make gdm syntax a little less arcane



## gdmize() takes an matrix or data.frame and adds a row for SampleID, which is
    # required by gdm. I usually have my rows named instead, hence the function.
    gdmize <- function(x, sampleids, sampleids_name="SampleID"){
        a <- data.frame(sampleids, as.data.frame(x))
        colnames(a)[1] <- sampleids_name
        return(a)
    }


## this is a function that makes a friendly gdm plot
    # x is a gdm model, result from function gdm()
    plot_gdm_jld <- function(x, points_color="darkslategray4", pred_colors="auto", 
        line_back_col="black", line_front_col="white", PSAMPLE=200, top_blank=FALSE,
        coef_threshold=0, points_cex=0.25, points_opacity=1){
        
        # for testing:
        # x <- mod_gdm; pred_color="auto"; coef_threshold=0; line_back_col="black"
        # line_front_col="white"; PSAMPLE=200; top_blank=FALSE; points_cex=1; points_color="darkslategray4"
        
        
        pc <- col2rgb(points_color) / 255
        points_color <- rgb(pc[1], pc[2], pc[3], alpha=points_opacity)



        require(gdm)
        # setting from original plot.gdm, not sure what it does so I'm leaving it
        options(warn.FPU = FALSE)

        ## define plot area type (2 cols, 1 row)
            par(mfrow=c(2,1), mai=c(1, 1, 0.1, 0.1))

        ## First plot - observed vs predicted compositional dissimilarity

            # make plot - blank if top_blank==T
            if(top_blank==TRUE){ptype<-"n"}else{ptype<-"p"}

            plot(x$predicted, x$observed, xlab = "Predicted community dissimilarity", 
                ylab = "Observed community dissimilarity",
                #ylim = c(0, 1),
                pch = 20, 
                cex = points_cex, 
                col = points_color,
                type = ptype
            )
            

            # add model fit
            if(top_blank==FALSE){
                overlayX <- overlayY <- seq(from = min(x$predicted), to = max(x$predicted), length = PSAMPLE)
                lines(overlayX, overlayY, lwd = 6, col=line_back_col)
                lines(overlayX, overlayY, lwd = 2, col=line_front_col, lty=2)
            }

        ## Organize spline data

            # figure out how many predictors we need to plot
            n_preds <- length(x$predictors)

            # make data frame for spline info (this makes code for extracting plot data 50000% more legible)
            spline_df <- data.frame(
                pred_ind=rep(1:length(x$predictors), x$splines),
                pred_name=rep(x$predictors, x$splines),
                coefficient=x$coefficients,
                knot=x$knots
            )

            # standardize all knots to 0-1 range
            # this allows all splines to be plotted together
            for(p in x$predictors){
                knots_i <- spline_df$knot[spline_df$pred_name == p]
                knots_i <- (knots_i - min(knots_i)) / (max(knots_i) - min(knots_i))
                spline_df$knot[spline_df$pred_name == p] <- knots_i
            }


            # make a list of predictor plot data
            pred_plot_list <- list()
            # fill list up 
            for(i in 1:n_preds){
                
                # not sure why this is pre-allocated, but it can't hurt much to leave it
                preddata_i <- rep(0, times = PSAMPLE)
                # c function to get predictor plot data. No idea what it does or how it works. 
                # I had to reverse-engineer the arguments it takes, but it works 100% now.
                pred_plot_list[[i]] <- .C("GetPredictorPlotData", 
                    pdata = as.double(preddata_i), 
                    as.integer(PSAMPLE), 
                    as.double(spline_df$coefficient[spline_df$pred_ind == i]),
                    as.double(spline_df$knot[spline_df$pred_ind == i]),
                    as.integer( sum(spline_df$pred_ind == i) ),
                    PACKAGE = "gdm"
                )
                # named lists are nice, add name 
                names(pred_plot_list)[i] <- x$predictors[i]
            }

            maxcoef <- sapply(X=pred_plot_list, FUN=function(x){max(x$pdata)})
            mmax <- max(maxcoef)
            # re-scale everything if something is > 1!
            if(mmax > 1){
                for(i in 1:n_preds){
                    pred_plot_list[[i]]$pdata <- pred_plot_list[[i]]$pdata / mmax
                    pred_plot_list[[i]][[3]] <- pred_plot_list[[i]][[3]] / mmax
                }
                maxcoef <- maxcoef / mmax
            }

            # drop variables that have max coef below threshold
            goodpreds <- x$predictors[maxcoef >= coef_threshold]
            pred_plot_list <- pred_plot_list[names(pred_plot_list) %in% goodpreds]
            n_preds <- length(pred_plot_list)

        ## Second plot - plot splines
            # make empty plot frame
            plot(x=NULL, y=NULL, type="n", 
                xlim=c(0, 1),
                ylim=c(0, max(maxcoef)),
                xlab="Variable dissimilarity",
                ylab="Partial community distance"
            )

            # make colors
            # if auto (default), make some colors
            if(pred_colors[1] == "auto"){
                # R interpreter checks logic first, so this is OK
                pred_colors <- rainbow(n_preds, start=0, end=0.60)
            }else{
                # make sure user colors are long enough
                while(length(pred_colors) < n_preds){
                    pred_colors <- c(pred_colors, pred_colors)
                }
                # trim
                pred_colors <- pred_colors[1:n_preds]
            }

            # plot 'em
            for(i in 1:length(pred_plot_list)){
                points(
                    x=seq(from=0, to=1, length=PSAMPLE),
                    y=pred_plot_list[[i]]$pdata,
                    type="l", 
                    col=pred_colors[i],
                    lwd=6
                )
            }

            # I don't know why this was hre?
            #for(i in 1:length(maxcoef)){
            #    maxcoef[i] <- round(max(pred_plot_list[[i]]$pdata), 2)
            #}
            # order for maxcoef hi2low
            hi2low_order <- order(maxcoef, decreasing = TRUE)
            legend_labels <- paste0(sprintf("%1.3f", round(maxcoef,3)), " - ", names(pred_plot_list), sep="")

            # add legend
            legend(x=0, y=max(maxcoef), legend=legend_labels[hi2low_order], col=pred_colors[hi2low_order], 
                lty=1, lwd=6, bty = "n")
    }

## function to reset par()
    resetPar <- function() {
        dev.new()
        op <- par(no.readonly = TRUE)
        dev.off()
    }

## visualize pairwise relationships among variables in a nice and minimalist way
    # r's pairs() isn't good enough for me, and chart.Correlation from performanceAnalytics is too messy
    # df is just a data frame where each column is a numeric variable to plot
    plot_pairwise_corrs <- function(df, label_cex=1, point_cex=1, cor_cex=2, cor_red_lim=0.70, mthd="pearson"){
        n <- ncol(df)
        par(mfrow = c(n,n), oma = c(5,4,0,0), mar = c(0,0,0,0) )
        # make a matrix to figure out which type of plot to do at position i,j
        # lower tri = scatterplots, diag=names, upper tri = correlation coefficients
        typemat <- matrix("D", nrow=n, ncol=n)
        typemat[lower.tri(typemat)] <- "L"
        typemat[upper.tri(typemat)] <- "U"
        for(i in 1:n){for(j in 1:n){
            if(typemat[i,j] == "L"){
                # lower tri - do scaterplot
                plot(x=df[,j], y=df[,i], axes = FALSE, xlab="", ylab="", pch=20, cex=point_cex)
                box()
            }else if(typemat[i,j] == "D"){
                # diag - write variable name
                plot(1, type="n", xlim=c(-1, 1), ylim=c(-1, 1), axes = FALSE, xlab="", ylab="", pch=20)
                text(x=0, y=0, labels=colnames(df)[i], cex=label_cex, srt=-45)
                box()
            }else if(typemat[i,j] == "U"){
                # upper tri - nicely display correlation coefficient (r)
                cor_ij <- cor(df[,j], df[,i], use="complete.obs", method=mthd) 
                if(cor_ij > cor_red_lim || cor_ij < (-1 * cor_red_lim)){
                    col_ij <- "red"
                }else{
                    col_ij <- "black"
                }
                cor_ij <- sprintf("%.2f", round(cor_ij,2))

                plot(1, type="n", xlim=c(-1, 1), ylim=c(-1, 1), axes = FALSE, xlab="", ylab="", pch=20)
                text(x=0, y=0, labels=cor_ij, cex=cor_cex, col=col_ij)
                box()
            }
        }}
        resetPar()
    }


#     spt2 <- site_pair_from_list( responseMat=unifrac_beta_distmat, predList=predictor_list2 )
#     responseMat <- unifrac_beta_distmat
#     predList <- predictor_list2
#     preds2use <- NULL


## site_pair_from_list generates GDM's sitepair table from a list of objects
    # valid types in the list are : "numeric", "matrix", or "list"
    # see example above
    # a strength of this approach is that one can use only a subset of all the predictors
    # with the preds2use argument.
    site_pair_from_list <- function(responseMat, predList, preds2use=NULL){
        # if preds2use is specified, simplify predList accordingly
        if(!is.null(preds2use)){
            # drop unused items from predList
            predList <- predList[names(predList) %in% preds2use]
            # get predList into the same order as preds2use (only matters for metadata column vectors...)
            predList <- predList[order(match(names(predList), preds2use))]
        }
        # get classes
        predClasses <- lapply(X=predList, FUN=class)

        # make table - if ONLY MATRIX, make data with fake variable instead.
        if(sum(predClasses == "numeric") > 0){
            predDF <- data.frame(
                SampleID=rownames(responseMat),                     # siteColumn
                simplify2array(predList[predClasses %in% c("integer", "numeric")]) # data columns
            )
        }else{
            predDF <- data.frame(
                SampleID=rownames(responseMat),                     # siteColumn
                FakeData=rep(0, nrow(responseMat))                  # fake data column
            )
        }

        # check to make sure "matrix_n" isn't used as a name for any predictors,
        # since GDM uses that internally for some stuff and I replace it later
        # on with better names. 
        if(any(grepl(x=names(predList), pattern="matrix_[0-9]+"))){
            stop("Predictors can't be named \"matrix_n\".")
        }

        # check if geo is included (one and only one time!)
        # if so, add geo information to predDF
        # if not, add fake geo information (because formatsitepair() is dumb)
        if(sum(predClasses == "list") == 1){
            geoLat <- predList[[which(predClasses=="list")]]$Lat
            geoLon <- predList[[which(predClasses=="list")]]$Lon
        }else{
            geoLat <- rep(1, nrow(responseMat))
            geoLon <- rep(1, nrow(responseMat))
        }
        # add real or fake lat/longs to predDF
        predDF <- data.frame(
            predDF,
            Lat=geoLat,
            Lon=geoLon
        )

        # format distance matrices
        if(sum(predClasses == "matrix") > 0){
            matrixList <- predList[predClasses == "matrix"]
            matrixList <- lapply(X=matrixList, FUN=gdmize, sampleids=rownames(responseMat))
        }else{
            matrixList <- NULL
        }
        
        # make sitepair table
        spt <- formatsitepair(
            bioData=gdmize(responseMat, rownames(responseMat)), bioFormat=3,
            predData=predDF,
            XColumn="Lon", YColumn="Lat",
            distPreds=matrixList,
            siteColumn="SampleID"
        )

        # remove NAs
        spt <- na.omit(spt)

        # fix matrix names to be original names
        mat_names <- names(matrixList)
        if(length(mat_names) > 0){for(i in 1:length(mat_names)){
            gdmname_i <- paste0("matrix_", i)
            colnames(spt) <- sub(x=colnames(spt), pattern=gdmname_i, replacement=mat_names[i])
        }}

        # remove "FakeData" cols from spt
        spt <- spt[, !grepl(pattern="FakeData", x=colnames(spt))]

        return(spt)
    }

## forward_adonis
    # forward model selection for adonis
    # all RHS vars must be column vectors, within a matrix.
    # no interaction terms are considered.
    # LHS is a dist object, maybe a community data matrix would work, not tested.
    fwd_adonis <- function(lhs, rhs, ncores=4){
        require(parallel)
        vars_in_model <- NULL
        lhs_name <- deparse(substitute(lhs))
        Ps <- R2s <- matrix(data=NA, nrow=ncol(rhs), ncol=ncol(rhs), dimnames=list(colnames(rhs)))

        # this function takes names of variables and returns a formula.
        makefrmla <- function(v, y="lhs"){ as.formula(paste(y, "~", paste(v, collapse=" + "))) }
        
        # start progress bar
        pb <- txtProgressBar(min=0, max=sum(1:ncol(R2s)), style=3)
        n_completed <- 0

        # do model selection
        for(j in 1:ncol(R2s)){
            newvars <- colnames(rhs)[! colnames(rhs) %in% vars_in_model]
            # get list of aov tables for each potential new model
            aovs_newvars <- mclapply(
                X=newvars, 
                FUN=function(x){ as.data.frame(adonis(makefrmla(c(vars_in_model, x)), data=rhs)$aov.tab) },
                mc.cores=ncores
            )
            # calculate total R2 for each potential model
            total_r2s <- sapply(
                X=aovs_newvars, 
                FUN=function(x){ sum(x$R2[! rownames(x) %in% c("Residuals", "Total")]) }
            )
            # choose which term to add based on total R2 of model
            toadd <- which.max(total_r2s)
            # add term to vars_in_model, and put Pvals and R2s in output matrices
            vars_in_model <- c(vars_in_model, newvars[toadd])
            for(i in 1:nrow(R2s)){
                term <- rownames(R2s)[i]
                if(term %in% rownames(aovs_newvars[[toadd]])){
                    R2s[i,j] <- round(aovs_newvars[[toadd]]$R2[ rownames(aovs_newvars[[toadd]]) == term ], 3)
                    Ps[i,j] <- round(aovs_newvars[[toadd]]$"Pr(>F)"[ rownames(aovs_newvars[[toadd]]) == term ], 3)
                }
            }
            # update progress bar
            n_completed <- n_completed + length(newvars)
            setTxtProgressBar(pb, n_completed)
        }
        message("")
        # which is the last model that had only significant variables?
        allsig <- apply(X=Ps, MAR=2, FUN=function(x){ all(x[!is.na(x)] < 0.05) })
        lastallsig <- which.max(which(allsig))
        if(length(lastallsig) <= 0){
            formula_sig <- "No significant models."
        }else{
            formula_sig <- makefrmla(v=vars_in_model[1:lastallsig], y=lhs_name)
        }
        formula_all <- makefrmla(v=vars_in_model, y=lhs_name)
        # return relevant objects
        return(list(
            formula_all,
            formula_sig,
            Pvals=Ps,
            R2s=R2s
        ))

    }
    

