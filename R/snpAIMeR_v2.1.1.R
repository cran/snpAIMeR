#' @title Assess the Diagnostic Power of Genomic Marker Combinations
#'
#' @description
#' Population genetics package for optimizing diagnostic panels. User-selected
#' candidate markers are assessed individually and in combination for how
#' well they can predict the source population of known samples. Requires a
#' genotype file in STRUCTURE format.
#'
#' @param run_mode Modes are "interactive", "non-interactive", or "example"; mode
#' must be in quotes.
#' @param config_file Yaml file required for "non-interactive" mode; filename/path
#' must be in quotes.
#' @param verbose Default is TRUE.
#'
#' @return Cross-validation assignment rates for individual markers, marker
#' combinations, and panel sizes. Outputs three .csv and two .pdf files to a
#' user-specified directory.
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#'
#' @details
#' Yaml file format for "non-interactive" mode (do not include bullet points):
#' \itemize{
#'   \item min_range: <minimum panel size>
#'   \item max_range: <maximum panel size; we recommend no more than 15 markers>
#'   \item assignment_rate_threshold: <value from 0 to 1>
#'   \item cross_validation_replicates: <we recommend 100 minimum>
#'   \item working_directory: <path name in quotes>
#'   \item structure_file: <path name in quotes>
#'   \item number_of_individuals: <same as adegenet's "n.ind">
#'   \item number_of_loci: <same as adegenet's "n.loc">
#'   \item one_data_row_per_individual: <TRUE or FALSE>
#'   \item column_sample_IDs: <column number>
#'   \item column_population_assignments: <column number>
#'   \item column_other_info: <column number>
#'   \item row_markernames: <row number>
#'   \item no_genotype_character: <default is "-9">
#'   \item optional_population_info: <optional>
#'   \item genotype_character_separator: <optional>
#' }
#'
#' Minimizing run time: Because of the number of possible combinations, we
#' recommend testing no more than 15 markers. For example, testing 15 markers
#' in panel sizes of 1 to 15 (32,767 total combinations) with 1,000
#' cross-validation replicates on a system with 48 processor cores took about 5
#' hours and 20 GB RAM. Reducing the number of cross-validation replicates will
#' reduce run time, however, we recommend no less than 100 replicates.
#'
#' @examples
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   data(nancycats, package = "adegenet")
#'   snpAIMeR("example", verbose = TRUE)
#' }
#'
#' @seealso https://github.com/OksanaVe/snpAIMeR


snpAIMeR <- function(run_mode, config_file = NULL, verbose = TRUE) {
  run_mode <- tolower(gsub("[^A-Za-z0-9]", "", run_mode))

  # Set up parallelization
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if ((nzchar(chk) && chk == "TRUE") || run_mode == "example") {
    # CRAN testing limits package examples to 2 cores
    n.cores <- 2L
  } else {
    # Setup backend to use multiple processors
    if (verbose == TRUE) {
      cat(parallel::detectCores(), "cores detected\n")
    }
    n.cores <- parallel::detectCores() - 1
  }
  # Create the cluster
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  if (verbose == TRUE) {
    print(my.cluster)
  }
  suppressWarnings(doParallel::registerDoParallel(cl = my.cluster))

  # Data input
  if((run_mode == "noninteractive") & (is.null(config_file))) {
    stop("Error: non-interactive mode requires config file")
  } else if ((run_mode == "noninteractive") & (!is.null(config_file))) {
      config = yaml::yaml.load_file(config_file)
      # snpAIMeR settings
      withr::local_dir(config$working_directory)
      min_loc_num <- as.integer(config$min_range)
      max_loc_num <- as.integer(config$max_range)
      assignment_rate_threshold <- as.double(config$assignment_rate_threshold)
      cv_replicates <- as.integer(config$cross_validation_replicates)
      # adegenet::read.structure settings
      pop_file <- config$structure_file
      number_of_individuals <- as.integer(config$number_of_individuals)
      number_of_loci <- as.integer(config$number_of_loci)
      column_sample_IDs <- as.integer(config$column_sample_IDs)
      column_population_assignments <- as.integer(config$column_population_assignments)
      row_markernames <- as.integer(config$row_markernames)
      column_other_info <- as.integer(config$column_other_info)
      if (verbose == TRUE) {
        pop_data <- adegenet::read.structure(pop_file, n.ind = number_of_individuals, n.loc = number_of_loci, onerowperind = config$one_data_row_per_individual, col.lab = column_sample_IDs, col.pop = column_population_assignments, col.others = column_other_info, row.marknames = row_markernames, NA.char = config$no_genotype_character, pop = config$optional_population_info, sep = config$genotype_character_separator, ask = FALSE, quiet = FALSE)
      } else {
          pop_data <- adegenet::read.structure(pop_file, n.ind = number_of_individuals, n.loc = number_of_loci, onerowperind = config$one_data_row_per_individual, col.lab = column_sample_IDs, col.pop = column_population_assignments, col.others = column_other_info, row.marknames = row_markernames, NA.char = config$no_genotype_character, pop = config$optional_population_info, sep = config$genotype_character_separator, ask = FALSE, quiet = TRUE)
      }
      loci = adegenet::locNames(pop_data)
      if (verbose == TRUE) {
        cat("Data file contains", length(loci), "markers\n")
        cat("File contains the following group definitions:")
        print(table(pop_data@pop))
      }
  } else if (run_mode == "interactive") {
      working_directory <- readline(prompt = "Enter path to working directory: ")
      withr::local_dir(working_directory)
      str_fn <- readline(prompt = "Enter path to STRUCTURE file: ")
      if (verbose == TRUE) {
        pop_data <- adegenet::read.structure(str_fn, quiet = FALSE)
      } else {
        pop_data <- adegenet::read.structure(str_fn, quiet = TRUE)
      }
      loci = adegenet::locNames(pop_data)
      if (verbose == TRUE) {
        cat("Data file contains", length(loci), "markers\n")
        cat("File contains the following group definitions:")
        print(table(pop_data@pop))
      }
      min_loc_num <- as.integer(readline(prompt = "Minimum number of markers in combination: "))
      max_loc_num <- as.integer(readline(prompt = "Maximum number of markers in combination: "))
      assignment_rate_threshold <- as.double(readline(prompt = "Threshold value for the minimum rate of successful assignments: "))
  	  cv_replicates <- as.integer(readline(prompt = "Number of cross-validation replicates: "))
  } else if (run_mode == "example") {
      # pop_data variable is replaced with nancycats
      # nancycats dataset is reduced from 9 to 3 markers
      nancycats_3_markers <- nancycats[,loc=c("fca8", "fca23", "fca43")]
      loci = adegenet::locNames(nancycats_3_markers)
      cat("Data file contains", length(loci), "markers\n")
      cat("File contains the following group definitions:")
      print(table(nancycats@pop))

      min_loc_num <- as.integer(1)
      max_loc_num <- as.integer(3)
      assignment_rate_threshold <- as.double(0.9)
      cv_replicates <- as.integer(5)
  } else {
      stop('\nError: invalid run_mode\n')
  }

  if (run_mode != "example") {
	 df = data.frame(marker = character(), avg_success_rate = double())
   utils::write.table(df,"Panel_size_assign_rate.csv", row.names = FALSE, sep = ",")
   utils::write.table(df,"All_combinations_assign_rate.csv", row.names = FALSE, sep = ",")
   utils::write.table(df,"Above_threshold_assign_rate.csv", row.names = FALSE, sep = ",")
  }
	mn_rate = data.frame(markers = integer(), mean_rate = double())

	if(min_loc_num == 1) {
		range = c(min_loc_num:max_loc_num)
	} else if(min_loc_num > 1) {
		range = c(1, min_loc_num:max_loc_num)
	} else if(min_loc_num < 1) {
		print("Incorrect number of markers provided")
	}

	for(n in range) {
		rate_ls <- list()
		loc_comb = utils::combn(loci,n)
		all_loci = data.frame(markers = character(), avg_success_rate = double())
		good_loci = data.frame(markers = character(), avg_success_rate = double())

		for(i in 1:ncol(loc_comb)){
		  if (run_mode == "example") {
		    test = nancycats[,loc=loc_comb[,i]]
		  } else {
		    test = pop_data[,loc = loc_comb[,i]]
		  }

			results <- foreach::foreach(j = 1:cv_replicates, .combine = 'c') %dopar% {
				tryCatch({
					kept.id <- unlist(tapply(1:adegenet::nInd(test), adegenet::pop(test), function(e) sample(e, round(0.75*length(e),0), replace = FALSE)))
					x <- test[kept.id]
					x.sup <- test[-kept.id]
					dapc4 <- adegenet::dapc(x, n.pca = 50, n.da = 20)
					pred.sup <- adegenet::predict.dapc(dapc4, newdata = x.sup)
					a <- mean(as.character(pred.sup$assign) == as.character(adegenet::pop(x.sup)))
				}, error = function(e){
				return(0)
			})
		}
			rt = list(results)
			rate_ls <- append(rate_ls, rt)

			markers = toString(loc_comb[,i])
			mean_rate <- sum(results)/length(results)
			if (verbose) {
			  print(paste0("Panel size ", n, ", combination ", i))
			  cat(markers, "\n")
			  cat(mean_rate, "\n")
			}
			all_loci[nrow(all_loci)+1,] = c(markers, mean_rate)
			if(mean_rate >= assignment_rate_threshold) {
				good_loci[nrow(good_loci)+1,] = c(markers, mean_rate)
			}
		}
		# Make box plot of each marker's cross-validation replicates
		if(n == 1) {
		  #Make the assignment rate tidy
		  loc_comb_list <- as.list(loc_comb)
		  loc_comb_df <- as.data.frame(do.call(rbind, loc_comb_list))
		  names(loc_comb_df)[1] <- "marker_name"
		  rate_df <- as.data.frame(do.call(rbind, rate_ls))
		  rate_df <- cbind(loc_comb_df, rate_df)
		  rate_df <- tidyr::pivot_longer(rate_df, -marker_name, names_prefix = "V", names_to = "replicate", values_to = "rate")
		  # Group the data by marker
		  rate_df_group_marker <- rate_df %>% dplyr::group_by(rate_df$marker_name)
		  # Get mean of all markers/replicates
		  mean_results <- round(mean(results), 4)
		  # Make box plot
		  box_plot_title <- paste0("Marker average = ", mean_results)
		  each_marker_plot <- ggplot2::ggplot(rate_df_group_marker, ggplot2::aes(x = forcats::fct_reorder(rate_df_group_marker$marker_name, readr::parse_number(rate_df_group_marker$marker_name)), y = rate_df_group_marker$rate)) +
		    ggplot2::geom_boxplot(fill = "gray", outlier.size = 0.5) +
		    ggplot2::labs(title = box_plot_title, x = "Marker", y = "Average cross-validation assignment rate") +
		    ggplot2::theme_light() +
		    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		  if (run_mode != "example") {
		    suppressMessages(ggplot2::ggsave("Single_marker_assign_rate.pdf"))
		  }
		  if (verbose == TRUE || run_mode == "example") {
		    print(each_marker_plot)
		  }

	   }
    # Histogram of rate distribution. Only prints to screen.
		else if(n > 1 && verbose == TRUE) {graphics::hist(results, main = paste0("Example cross-validation for ", n, " markers, " , cv_replicates, " replicates"), xlab = "Cross-validation assignment rate", ylab = "Frequency")
		  graphics::abline(v = mean(results), col = "red", lwd = 2)
		  graphics::mtext(paste("Mean =", round(mean(results), 4)), side = 3, col = "red")
		}
		# Get the mean assignment rate for each size group
		mn_rate[nrow(mn_rate)+1,] = c(n, mean(results))
		combination_mean_matrix <- cbind(mn_rate$markers, mn_rate$mean_rate)
		colnames(combination_mean_matrix) <- c("panel_size","avg_success_rate")
		combination_mean_df <- as.data.frame(combination_mean_matrix)

		if (run_mode != "example") {
		 utils::write.table(combination_mean_df, "Panel_size_assign_rate.csv", row.names = FALSE, col.names = TRUE, append = FALSE, sep = ",")
		 utils::write.table(all_loci, "All_combinations_assign_rate.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
		 utils::write.table(good_loci, "Above_threshold_assign_rate.csv", row.names = FALSE, col.names = FALSE, append = TRUE, sep = ",")
		}
	}
  # Make line plot of combination assignment rate means
	if (nrow(combination_mean_df) > 1) {
	  combination_means <- ggplot2::ggplot(combination_mean_df, ggplot2::aes(x = panel_size, y = avg_success_rate)) +
	    ggplot2::geom_line() +
	    ggplot2::geom_point() +
	    ggplot2::scale_x_continuous(breaks = seq(round(max(combination_mean_df$panel_size),0))) +
	    ggplot2::labs(title = "Panel sizes", x = "Number of markers in combination", y = "Average assignment rate") +
	    ggplot2::theme_light()
	   if (run_mode != "example"){
	     suppressMessages(ggplot2::ggsave("Panel_size_assign_rate.pdf"))
	   }
	   if (verbose == TRUE) {
	     print(combination_means)
	   }
	}
	if (verbose == TRUE && run_mode != "example") {
	  cat(nrow(good_loci), "marker combinations passed the assignment rate threshold", assignment_rate_threshold, "\n")
	}

	suppressWarnings(parallel::stopCluster(cl = my.cluster))
}
