library("ReportingTools")	# For creating HTML reports
library("knitr")		# For creating HTML reports

library("lattice") # For plotting


# ================ Hard coded HTML code changes =================
add_index_col_fct <- "
  [...document.querySelectorAll('#id2 tr')].forEach((row, i) => {
    var cell = document.createElement(i<2 ? 'th' : 'td');

    if (i <2) {
        row.insertCell(0);
    } else {
        var cell = row.insertCell(0);
        cell.classList.add('my_index_col');
        cell.innerHTML = (i-1);
    }
});"
add_index_col_update_fct <- "
  t.on('draw.dt', function(){
    console.log('Update index');
    let n = 0;
    $('.my_index_col').each(function () {
        $(this).html(++n);
    })
})"


# 1. Insert index column -> Must be inserted before DataTable is initialized
original_table_init_fct_head <- "function configureTable(i, el) {"
substitute_table_init_fct_head <- paste(original_table_init_fct_head, add_index_col_fct)

# 2. Remove pre-ordering of table
remove_to_disable_preordering <- "\"aaSorting\":[[0,'asc']],"

# ------------- Following hack code is not needed anymore... -----------------
# 3. Create local variable "t" that references the datatable
original_js_fct_head <- "$(this).dataTable({"
substitute_js_fct_head <- paste("var t = ", original_js_fct_head)

# 4. Use this as anchor to add more JS code
original_js_fct_tail <- '}).columnFilter({sPlaceHolder: "head:before",
                                aoColumns : filterClasses
                                });'



create_html_table <- function(input_file_path, sep="\t", title="Report Title", info_text="Info Text",
                              base_name="my_report", output_dir=".") {
  # Load table from data file
  # as.data.frame(resOrdered)
  input_table <- as.data.frame(read.csv(input_file_path, header=TRUE, sep=sep))
  if (nrow(input_table) == 0) {
    input_table[nrow(input_table)+1,] <- "No data"
  }
  # Remove column with no header (R names them "X")
  remove.cols <- names(input_table) %in% c("", "X")
  input_table <- input_table[! remove.cols]

  # Use ReportingTools to automatically generate dynamic HTML documents
  html_report <- ReportingTools::HTMLReport(shortName=base_name, title=title,
                                            reportDirectory=output_dir)

  # 1. Add a table to the report
  ReportingTools::publish(input_table, html_report)
  # 2. Add info text to the report
  ReportingTools::publish(info_text, html_report)

  # Also graphs can be added to the report
  # # Randomly
  # y <- rnorm(500)
  # plot<-lattice::histogram(y, main="Sample of 500 observations from a Normal (0,1)")
  # # 3. Add plot to the report
  # ReportingTools::publish(plot, html_report)

  # Finally, create the report
  ReportingTools::finish(html_report)
}

replace_external_scripts_and_styles <- function(input_html_file) {
  # Replace external scripts and styles with internal copies

  # Read input file
  html_file_content <- readLines(input_html_file, warn=FALSE)

  external_js_scripts <- c('<script language="JavaScript" src="jslib/jquery-1.8.0.min.js"></script>',
    '<script language="JavaScript" src="jslib/jquery.dataTables-1.9.3.js"></script>',
    '<script language="JavaScript" src="jslib/bootstrap.js"></script>',
    '<script language="JavaScript" src="jslib/jquery.dataTables.columnFilter.js"></script>',
    '<script language="JavaScript" src="jslib/jquery.dataTables.plugins.js"></script>',
    '<script language="JavaScript" src="jslib/jquery.dataTables.reprise.js"></script>',
    '<script language="JavaScript" src="jslib/bootstrap.js"></script>')

  external_css_styles <- c('<link rel="stylesheet" type="text/css" href="csslib/bootstrap.css" />',
    '<link rel="stylesheet" type="text/css" href="csslib/reprise.table.bootstrap.css" />')


  # Replace external scripts with CDN versions
  jquery_cdn <- '<script src="https://code.jquery.com/jquery-3.6.3.min.js" integrity="sha256-pvPw+upLPUjgMXY0G+8O0xUf+/Im1MZjXxxgOcBQBXU=" crossorigin="anonymous"></script>'
  jquery_datatable_cdn <- '<script src="https://cdn.datatables.net/1.13.1/js/jquery.dataTables.min.js"></script>'
  bootstrap_js_cdn <- '<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/js/bootstrap.min.js"></script>'
  html_file_content <- sub(external_js_scripts[1], jquery_cdn, html_file_content)
  html_file_content <- sub(external_js_scripts[2], jquery_datatable_cdn, html_file_content)
  html_file_content <- sub(external_js_scripts[3], bootstrap_js_cdn, html_file_content)

  # Replace external styles with CDN versions
  bootstrap_css_cdn <- '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css" />'
  html_file_content <- sub(external_css_styles[1], bootstrap_css_cdn, html_file_content)

  # Replace external scripts with local copies
  for (js_import in external_js_scripts[4:length(external_js_scripts)]) {
    js_source_file <- sub('<script language="JavaScript" src="', '', js_import)
    js_source_file <- sub('"></script>', '', js_source_file)
    path_to_js_file <- file.path(dirname(input_html_file), js_source_file)
    js_code <- paste(readLines(path_to_js_file, warn=FALSE), collapse="\n")

    # Replace external script with internal script
    html_file_content <- sub(js_import, paste('<script language="JavaScript">', js_code, '</script>', sep="\n"), html_file_content)
  }

  # Replace external styles with local copies
  for (css_import in external_css_styles[2:length(external_css_styles)]) {
    css_source_file <- sub('<link rel="stylesheet" type="text/css" href="', '', css_import)
    css_source_file <- sub('" />', '', css_source_file)
    path_to_css_file <- file.path(dirname(input_html_file), css_source_file)
    css_code <- paste(readLines(path_to_css_file, warn=FALSE), collapse="\n")

    # Replace external style with internal style
    html_file_content <- sub(css_import, paste('<style>', css_code, '</style>', sep="\n"), html_file_content)
  }

  return(html_file_content)
  # Write HTML file
}

add_index_column_functionality <- function(input_html_content) {
  "
  Add index column functionality to the HTML table
  Also disable pre-sorting by the first column
  "
  # 1. Insert index column -> Must be inserted before DataTable is initialized
  input_html_content <- sub(original_table_init_fct_head, substitute_table_init_fct_head, input_html_content, fixed=TRUE)
  # 2. Remove pre-ordering of table
  input_html_content <- sub(remove_to_disable_preordering, "", input_html_content, fixed=TRUE)

  return(input_html_content)
}

add_csv_download_button <- function(input_html_content) {
  "
  Add CSV download button to the HTML table.
  Button has class 'buttons-csv'.
  "
  # DataTables: Select only CSV button in intialization
  original_initialization <- "$(this).dataTable({"
  new_initialization <- "$(this).dataTable({\n\"buttons\": [\"csvHtml5\"],"

  # DataTables: Integration of buttons into DOM
  original_dom_declaration <- "\"sDom\": \"<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>\","
  new_dom_declaration <- "\"sDom\": \"<'row'<'span6'lB><'span6'f>r>t<'row'<'span6'i><'span6'p>>\","

  # JS libraries
  additional_js_lib_1 <- '<script src="https://cdn.datatables.net/buttons/2.3.6/js/dataTables.buttons.min.js"></script>'
  additional_js_lib_2 <- '<script src="https://cdn.datatables.net/buttons/2.3.6/js/buttons.html5.min.js"></script>'

  # CSS changes
  # Make position relative and float right
  additional_css_changes <- "<style> .buttons-csv { position: relative; float: right; } </style>"

  # Button classes
  # -> Add btn-primary class to CSV button (which has class 'buttons-csv')
  add_class_script <- "<script>$(document).ready(function(){$('button.buttons-csv').addClass('btn btn-sm btn-primary mb-2');} );</script>"

  # Introduce changes
  # 0. DataTables initialization
  input_html_content <- sub(original_initialization, new_initialization, input_html_content, fixed=TRUE)
  # 1. DOM declaration
  input_html_content <- sub(original_dom_declaration, new_dom_declaration, input_html_content, fixed=TRUE)
  # 2. JS libraries
  input_html_content <- sub("</head>", paste(additional_js_lib_1, additional_js_lib_2, "</head>", sep="\n"), input_html_content, fixed=TRUE)
  # 3. CSS changes
  input_html_content <- sub("</head>", paste(additional_css_changes, "</head>", sep="\n"), input_html_content, fixed=TRUE)
  # 4. Button classes
  input_html_content <- sub("</body>", paste(add_class_script, "</body>", sep="\n"), input_html_content, fixed=TRUE)

  return(input_html_content)
}

fix_table_width <- function(input_html_content) {
  "
    Fix table width to 100% -> Make it scrollable
  "
  # Insert wrapper at initialization to manage scrolling (scrollX has issue with alignment of headers...)
  original_initialization <- "$(this).dataTable({"
  new_initialization <- paste(original_initialization, '"initComplete": function (settings, json) {
      $(this).wrap("<div style=\'overflow:auto; width:100%; position:relative;\'></div>");
    },', sep="\n")

  # Ellipsis style for long text
  additional_css_changes <- "<style> table.dataTable td  {
        max-width: 250px;
        white-space: nowrap;
        text-overflow: ellipsis;
        overflow: hidden;
      }
    </style>"

  # 1. DataTables initialization
  input_html_content <- sub(original_initialization, new_initialization, input_html_content, fixed=TRUE)
  # 2. CSS changes
  input_html_content <- sub("</head>", paste(additional_css_changes, "</head>", sep="\n"), input_html_content, fixed=TRUE)

  return(input_html_content)
}


# Main function
main <- function() {
  # Import snakemake arguments
  input_files <- snakemake@params[["input_files"]]
  input_separators <- snakemake@params[["data_separators"]]
  input_titles <- snakemake@params[["data_titles"]]
  input_info_texts <- snakemake@params[["info_texts"]]
  output_dir <- snakemake@params[["html_output_dir"]]
  output_file_basenames <- snakemake@params[["html_output_file_basenames"]]

  # Iterate over all inputs and create HTML-reports
  for (i in 1:length(input_files)) {
    output_file_basename <- output_file_basenames[i]
    output_html_file <- file.path(output_dir, output_file_basename)

    print(paste("Creating report for", input_files[i]))
    print(paste("Output file basename:", output_file_basename))
    create_html_table(input_files[i], sep=input_separators[i], title=input_titles[i],
                      info_text=input_info_texts[i], base_name=output_file_basename,
                      output_dir=output_dir)

    # Replace external scripts and styles with internal copies
    updated_html_file <- replace_external_scripts_and_styles(output_html_file)

    # Add index column functionality
    updated_html_file <- add_index_column_functionality(updated_html_file)

    # Add CSV download button
    updated_html_file <- add_csv_download_button(updated_html_file)

    # Fix table width
    updated_html_file <- fix_table_width(updated_html_file)

    # Write updated HTML file
    writeLines(updated_html_file, output_html_file)
  }
}

main()
