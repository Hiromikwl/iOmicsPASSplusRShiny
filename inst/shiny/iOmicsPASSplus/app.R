# iOmicsPASS+ Shiny GUI
# Requires: shiny, bslib, DT, iOmicsPASSplus

library(shiny)
library(bslib)
library(DT)
library(visNetwork)
library(iOmicsPASSplus)

options(shiny.maxRequestSize = 2 * 1024^3)  # 2 GB upload limit

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ── HGNC reference gene lists (bundled) ───────────────────────────────────────
.hgnc_all <- {
  f <- system.file("extdata", "hgnc_genes.rds", package = "iOmicsPASSplus")
  if (nzchar(f) && file.exists(f)) unique(readRDS(f)$symbol) else character(0)
}
.hgnc_protein_product <- {
  f <- system.file("extdata", "hgnc_protein_product.rds", package = "iOmicsPASSplus")
  if (nzchar(f) && file.exists(f)) readRDS(f) else character(0)
}

# ── Session working directory ──────────────────────────────────────────────────
# Created once per R process; survive across multiple browser connections to the
# same server instance (acceptable for a desktop app launched with run_iOmicsPASS()).
app_dir <- file.path(tempdir(), "iOmicsPASS_shiny_session")
dir.create(file.path(app_dir, "iOmicsPASS", "inputFiles"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_dir, "iOmicsPASS", "Output"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_dir, "Output"),                    recursive = TRUE, showWarnings = FALSE)

# ── Helpers ────────────────────────────────────────────────────────────────────
read_upload <- function(path) {
  df <- tryCatch(
    read.delim(path, sep = "\t", as.is = TRUE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) <= 1) {
    df <- tryCatch(
      read.csv(path, as.is = TRUE, check.names = FALSE),
      error = function(e) NULL
    )
  }
  df
}

ok_badge <- function(msg) {
  tags$span(style = "color:#198754; font-weight:500;",
            "✔ ", msg)
}
warn_badge <- function(msg) {
  tags$span(style = "color:#dc3545; font-weight:500;",
            "⚠ ", msg)
}
info_badge <- function(msg) {
  tags$span(style = "color:#6c757d;", "ℹ ", msg)
}

# ── Theme ──────────────────────────────────────────────────────────────────────
app_theme <- bs_theme(
  version   = 5,
  bg        = "#f8f9fa",
  fg        = "#212529",
  primary   = "#1a6faf",
  success   = "#198754",
  warning   = "#ffc107",
  danger    = "#dc3545",
  "navbar-bg" = "#1a6faf",
  "navbar-light-color" = "#ffffff",
  "navbar-light-active-color" = "#ffffff",
  "navbar-nav-link-font-size" = "1.05rem",
  "navbar-nav-link-padding-x" = "1rem"
)

# ── UI ─────────────────────────────────────────────────────────────────────────
.reset_file_js <- tags$script(HTML('
  Shiny.addCustomMessageHandler("resetFileInput", function(id) {
    var el = document.getElementById(id);
    if (el && el.type === "file") { el.value = ""; $(el).trigger("change"); }
    var prog = document.getElementById(id + "_progress");
    if (prog) { $(prog).hide().find(".progress-bar").css("width", "0"); }
  });

  Shiny.addCustomMessageHandler("exportVisNetwork", function(msg) {
    var container = document.getElementById(msg.id);
    if (!container) { console.warn("exportVisNetwork: container not found", msg.id); return; }
    var SCALE = 3;

    function doDownload(dUrl) {
      var link = document.createElement("a");
      link.href = dUrl; link.download = msg.name + "." + msg.type;
      document.body.appendChild(link); link.click(); document.body.removeChild(link);
      if (Shiny && Shiny.setInputValue)
        Shiny.setInputValue("net_canvas_b64", dUrl, {priority: "event"});
    }

    // Try to get vis-network + legend instances
    var network = null, legendNetwork = null;
    try {
      var wi = HTMLWidgets.getInstance(container);
      if (wi) {
        network       = wi.network || (wi.instance && wi.instance.network) || null;
        legendNetwork = wi.legend  || (wi.instance && wi.instance.legend)  || null;
      }
    } catch(e) {}
    if (!legendNetwork) {
      var legEl = document.getElementById("legend" + msg.id);
      if (legEl) legendNetwork = legEl.network || null;
    }

    // Snapshot canvas positions BEFORE any setSize distorts getBoundingClientRect
    var allCvs = Array.prototype.slice.call(container.querySelectorAll("canvas"))
                      .filter(function(c){ return c.width > 0 && c.height > 0; });
    if (!allCvs.length) { console.warn("exportVisNetwork: no canvases"); return; }
    var canvasInfo = allCvs.map(function(c) {
      var r = c.getBoundingClientRect();
      return { canvas: c, origW: r.width, origH: r.height };
    });
    var mainInfo = canvasInfo.reduce(function(best, ci) {
      return ci.canvas.width * ci.canvas.height > best.canvas.width * best.canvas.height ? ci : best;
    }, canvasInfo[0]);

    if (network && network.setSize) {
      var origW = container.offsetWidth;
      var origH = container.offsetHeight;
      var expW  = origW * SCALE;
      var expH  = origH * SCALE;

      network.setSize(expW + "px", expH + "px");
      // Also resize the legend vis-network so it re-renders at 3× resolution
      if (legendNetwork && legendNetwork.setSize) {
        canvasInfo.forEach(function(ci) {
          if (ci !== mainInfo)
            legendNetwork.setSize((ci.origW * SCALE) + "px", (ci.origH * SCALE) + "px");
        });
      }

      setTimeout(function() {
        network.redraw();
        if (legendNetwork && legendNetwork.redraw) legendNetwork.redraw();
        setTimeout(function() {
          var out = document.createElement("canvas");
          out.width = expW; out.height = expH;
          var ctx = out.getContext("2d");
          ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high";
          ctx.fillStyle = msg.background || "white"; ctx.fillRect(0, 0, expW, expH);
          // Main canvas fills the full output
          ctx.drawImage(mainInfo.canvas, 0, 0, expW, expH);
          // Legend canvases: pin top-right at 3× their original CSS size
          canvasInfo.forEach(function(ci) {
            if (ci === mainInfo) return;
            var lW = Math.round(ci.origW * SCALE);
            var lH = Math.round(ci.origH * SCALE);
            ctx.drawImage(ci.canvas, expW - lW, 0, lW, lH);
          });
          var mime = (msg.type === "jpeg") ? "image/jpeg" : "image/png";
          doDownload(out.toDataURL(mime, 0.95));
          // Restore both networks to original size
          network.setSize(origW + "px", origH + "px");
          if (legendNetwork && legendNetwork.setSize) {
            canvasInfo.forEach(function(ci) {
              if (ci !== mainInfo)
                legendNetwork.setSize(ci.origW + "px", ci.origH + "px");
            });
          }
          setTimeout(function() {
            network.redraw();
            if (legendNetwork && legendNetwork.redraw) legendNetwork.redraw();
          }, 80);
        }, 350);
      }, 80);
    } else {
      // Fallback: no vis-network instance; capture at 2× DPR scale
      var dpr  = Math.max(window.devicePixelRatio || 1, 2);
      var outW = (container.offsetWidth  || 800) * dpr;
      var outH = (container.offsetHeight || 600) * dpr;
      var out = document.createElement("canvas");
      out.width = outW; out.height = outH;
      var ctx = out.getContext("2d");
      ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high";
      ctx.fillStyle = msg.background || "white"; ctx.fillRect(0, 0, outW, outH);
      ctx.drawImage(mainInfo.canvas, 0, 0, outW, outH);
      canvasInfo.forEach(function(ci) {
        if (ci === mainInfo) return;
        var lW = Math.round(ci.origW * dpr), lH = Math.round(ci.origH * dpr);
        ctx.drawImage(ci.canvas, outW - lW, 0, lW, lH);
      });
      var mime = (msg.type === "jpeg") ? "image/jpeg" : "image/png";
      doDownload(out.toDataURL(mime, 0.95));
    }
  });
'))

ui <- page_navbar(
  title = tags$a(
    href = "#",
    style = "display:flex; align-items:center; text-decoration:none; height:100%;",
    tags$img(src = "logo.png", height = "115px",
             style = "display:block; filter:drop-shadow(0 1px 2px rgba(0,0,0,.35));")
  ),
  theme      = app_theme,
  window_title = "iOmicsPASS+",
  bg         = "#1a6faf",
  header     = tagList(
    .reset_file_js,
    tags$style(HTML("
      #interactive_net .vis-legend {
        top: 10px !important;
        bottom: auto !important;
      }
    "))
  ),

  # ──────────────────────────────────────────────────────────────────────────
  # Tab 1: Welcome
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "Welcome",
    icon = icon("home"),
    # Center the navbar brand/logo vertically and style nav links larger.
    tags$head(tags$style(HTML("
      .navbar-brand { display:flex !important; align-items:center; padding-top:0; padding-bottom:0; }
      .navbar-nav .nav-link { font-size:1.05rem !important; padding-left:1rem !important; padding-right:1rem !important; }
      .navbar-nav .nav-link .fa, .navbar-nav .nav-link svg { font-size:1.1em; }
    "))),
    # Patch DT 0.34 bug: data.lazyRender throws when data.x===null on a hidden element.
    # tags$head() is hoisted to the page <head> by htmltools regardless of nesting depth.
    tags$head(tags$script(HTML("
(function() {
  function patchDT() {
    if (!window.Shiny || !window.Shiny.outputBindings) return false;
    var dtEntry = window.Shiny.outputBindings.bindingNames['datatables'];
    if (!dtEntry || !dtEntry.binding || !dtEntry.binding.renderValue) return false;
    if (dtEntry.binding.__dtNullPatched) return true;
    var _origAdapter = dtEntry.binding.renderValue;
    dtEntry.binding.renderValue = function(el, data) {
      if (!data || data.x === null || data.x === undefined) {
        if (window.jQuery) jQuery(el).empty().append('&nbsp;');
        return;
      }
      return _origAdapter.call(this, el, data);
    };
    dtEntry.binding.__dtNullPatched = true;
    return true;
  }
  if (!patchDT()) {
    document.addEventListener('DOMContentLoaded', function() {
      if (!patchDT()) setTimeout(patchDT, 50);
    });
  }
})();
"))),
    div(
      class = "py-4",
      style = "margin-left:208px; max-width:860px;",
      h2("Welcome to iOmicsPASS+"),
      p(class = "lead",
        "iOmicsPASS+ integrates up to three omics datasets using known biological networks ",
        "(TF-regulatory network, PPI networks, miRNA–gene target network, etc) or network inference methods ",
        "to identify subnetwork signatures that discriminate between phenotypic groups of interest. ",
        "Please read ", tags$a("Koh et al.", href = "https://www.nature.com/articles/s41540-019-0099-y",
          target = "_blank"), " for more information."),
      hr(),
      h4("Workflow Overview"),
      div(style = "overflow-x:auto; padding-bottom:4px;",
        div(
          style = "display:flex; flex-wrap:nowrap; align-items:stretch; justify-content:center; gap:6px; min-width:fit-content; margin:1rem auto;",
          div(class = "card text-center text-white shadow-sm",
              style = "width:130px; flex-shrink:0; height:175px; display:flex; flex-direction:column; align-items:center; justify-content:center; padding:10px 8px; background:linear-gradient(135deg,#1a6faf,#2196f3); border:none; border-radius:12px;",
              tags$div(style = "font-size:2.2rem; margin-bottom:6px; color:#ffe082;", icon("upload")),
              tags$b("1. Data Upload"),
              tags$div(style="opacity:.9; font-size:.73rem; margin-top:4px;", "Up to 3 types of omics data + Phenotype")),
          div(style = "font-size:1.6rem; color:#888; font-weight:700; flex-shrink:0; display:flex; align-items:center; padding:0 2px;", "→"),
          div(class = "card text-center text-white shadow-sm",
              style = "width:130px; flex-shrink:0; height:175px; display:flex; flex-direction:column; align-items:center; justify-content:center; padding:10px 8px; background:linear-gradient(135deg,#0d7377,#14a085); border:none; border-radius:12px;",
              tags$div(style = "font-size:2.2rem; margin-bottom:6px; color:#b2ebf2;", icon("sitemap")),
              tags$b("2. Network Setup"),
              tags$div(style="opacity:.9; font-size:.73rem; margin-top:4px;", "Known biological networks or glasso network inference")),
          div(style = "font-size:1.6rem; color:#888; font-weight:700; flex-shrink:0; display:flex; align-items:center; padding:0 2px;", "→"),
          div(class = "card text-center text-white shadow-sm",
              style = "width:130px; flex-shrink:0; height:175px; display:flex; flex-direction:column; align-items:center; justify-content:center; padding:10px 8px; background:linear-gradient(135deg,#117a65,#27ae60); border:none; border-radius:12px;",
              tags$div(style = "font-size:2.2rem; margin-bottom:6px; color:#a5d6a7;", icon("cogs")),
              tags$b("3. iOmicsPASS"),
              tags$div(style="opacity:.9; font-size:.73rem; margin-top:4px;", "Predictive subnetwork analysis")),
          div(style = "font-size:1.6rem; color:#888; font-weight:700; flex-shrink:0; display:flex; align-items:center; padding:0 2px;", "→"),
          div(class = "card text-center text-white shadow-sm",
              style = "width:130px; flex-shrink:0; height:175px; display:flex; flex-direction:column; align-items:center; justify-content:center; padding:10px 8px; background:linear-gradient(135deg,#e65100,#f57c00); border:none; border-radius:12px;",
              tags$div(style = "font-size:2.2rem; margin-bottom:6px; color:#ffe0b2;", icon("chart-bar")),
              tags$b("4. Results"),
              tags$div(style="opacity:.9; font-size:.73rem; margin-top:4px;", "Predictive network edges and Sample classification")),
          div(style = "font-size:1.6rem; color:#888; font-weight:700; flex-shrink:0; display:flex; align-items:center; padding:0 2px;", "→"),
          div(class = "card text-center text-white shadow-sm",
              style = "width:130px; flex-shrink:0; height:175px; display:flex; flex-direction:column; align-items:center; justify-content:center; padding:10px 8px; background:linear-gradient(135deg,#6c3483,#9b59b6); border:none; border-radius:12px;",
              tags$div(style = "font-size:2.2rem; margin-bottom:6px; color:#e1bee7;", icon("flask")),
              tags$b("5. Prediction"),
              tags$div(style="opacity:.9; font-size:.73rem; margin-top:4px;", "Classification of new dataset"))
        )
      ),
      hr(),
      h4("Quick Start"),
      p("Load the example bundled datasets (Tulip omics datasets: Proteomics, microRNA; Network files: ",
        "PPI, miRNA-gene target networks, Human Pathway annotations) to explore the full workflow. ",
        "More details on Tulip study can be found on ",
        tags$a("Choi et al.", href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6460474/",
          target = "_blank"), "."),
      actionButton("load_example", "Load Example Data",
                   class = "btn btn-success btn-lg",
                   icon = icon("database")),
      uiOutput("example_loaded_msg"),
      hr(),
      h4(icon("code-branch"), " R Package"),
      p("iOmicsPASS+ is also available as an R package hosted on GitHub. ",
        "You can install it directly in R using ",
        tags$code('devtools::install_github("cssblab/iOmicsPASSplus")'), "."),
      tags$a(
        href = "https://github.com/cssblab/iOmicsPASSplus",
        target = "_blank",
        class = "btn btn-outline-dark btn-sm",
        icon("github"), " cssblab/iOmicsPASSplus"
      )
    )
  ),

  # ──────────────────────────────────────────────────────────────────────────
  # Tab 2: Data Upload
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "Data Upload",
    icon = icon("upload"),
    layout_sidebar(
      sidebar = sidebar(
        width = 270,
        open  = TRUE,
        h5("Format Notes"),
        tags$ul(
          tags$li(tags$b("Omics files:"), " features as rows, samples as columns; first column = feature name; data values start from the second column."),
          tags$li(tags$b("Phenotype:"), " column 1 = sample labels, column 2 = Group, columns 3+ = covariates for class priors."),
          tags$li("Tab-delimited (.txt) or comma-delimited (.csv) accepted.")
        ),
        hr(),
        h5("Data Options"),
        checkboxInput("log_X",   "Log-transform Data X",       value = FALSE),
        checkboxInput("log_Y",   "Log-transform Data Y",       value = FALSE),
        checkboxInput("log_Z",   "Log-transform Data Z",       value = FALSE),
        checkboxInput("std_data","Standardize data (recommended)", value = TRUE)
      ),

      div(
        class = "container-fluid py-3",
        # ── Data X ───────────────────────────────────────────────────────────
        card(
          class = "mb-3",
          card_header(tags$b("Data X"), tags$small(" (required)", style = "color:gray;")),
          div(class = "row g-3",
            div(class = "col-md-4",
              p(tags$small("Features in rows, samples in columns. Data values start from the second column.")),
              fileInput("file_X", NULL, placeholder = "Upload Data X (.txt/.csv)",
                        accept = c(".txt",".csv")),
              div(class = "d-flex justify-content-end mb-1",
                  actionButton("reset_X", "Reset", icon = icon("times"),
                               class = "btn btn-sm btn-outline-secondary")),
              uiOutput("status_X")
            ),
            div(class = "col-md-8",
              strong("Preview:"),
              div(style = "font-size:0.85em;",
                  DT::dataTableOutput("preview_X_tbl"))
            )
          )
        ),
        # ── Data Y ───────────────────────────────────────────────────────────
        card(
          class = "mb-3",
          card_header(tags$b("Data Y"), tags$small(" (optional)", style = "color:gray;")),
          checkboxInput("use_Y", "Include Data Y", value = FALSE),
          conditionalPanel(
            condition = "input.use_Y",
            div(class = "row g-3",
              div(class = "col-md-4",
                p(tags$small("Features in rows, samples in columns. Data values start from the second column.")),
                fileInput("file_Y", NULL, placeholder = "Upload Data Y (.txt/.csv)",
                          accept = c(".txt",".csv")),
                div(class = "d-flex justify-content-end mb-1",
                    actionButton("reset_Y", "Reset", icon = icon("times"),
                                 class = "btn btn-sm btn-outline-secondary")),
                uiOutput("status_Y")
              ),
              div(class = "col-md-8",
                strong("Preview:"),
                div(style = "font-size:0.85em;",
                    DT::dataTableOutput("preview_Y_tbl"))
              )
            )
          )
        ),
        # ── Data Z ───────────────────────────────────────────────────────────
        card(
          class = "mb-3",
          card_header(tags$b("Data Z"), tags$small(" (optional)", style = "color:gray;")),
          checkboxInput("use_Z", "Include Data Z", value = FALSE),
          conditionalPanel(
            condition = "input.use_Z",
            div(class = "row g-3",
              div(class = "col-md-4",
                p(tags$small(
                  "Features in rows, samples in columns. Data values start from the second column.",
                  tags$br(),
                  tags$b("Note:"), " Used as a normalization factor."
                )),
                fileInput("file_Z", NULL, placeholder = "Upload Data Z (.txt/.csv)",
                          accept = c(".txt",".csv")),
                radioButtons("normalizeZby", "Normalize Data Z by:",
                             choices = c("X" = "X", "Y" = "Y"), inline = TRUE),
                div(class = "d-flex justify-content-end mb-1",
                    actionButton("reset_Z", "Reset", icon = icon("times"),
                                 class = "btn btn-sm btn-outline-secondary")),
                uiOutput("status_Z")
              ),
              div(class = "col-md-8",
                strong("Preview:"),
                div(style = "font-size:0.85em;",
                    DT::dataTableOutput("preview_Z_tbl"))
              )
            )
          )
        ),
        # ── Row 2: Phenotype File, Pathway Annotation ─────────────────────────
        fluidRow(
          column(6,
            card(
              card_header(tags$b("Phenotype File"), tags$small(" (required)", style = "color:gray;")),
              p(tags$small(
                "Requires at least 2 columns. ",
                tags$b("Column 1:"), " sample labels. ",
                tags$b("Column 2:"), " Group. ",
                tags$b("Columns 3+:"), " covariates to be included for calculating class priors."
              )),
              fileInput("file_pheno", NULL, placeholder = "Upload Phenotype (.txt/.csv)",
                        accept = c(".txt",".csv")),
              div(class = "d-flex justify-content-end mb-1",
                  actionButton("reset_pheno", "Reset", icon = icon("times"),
                               class = "btn btn-sm btn-outline-secondary")),
              uiOutput("status_pheno"),
              tags$hr(),
              strong("Preview:"),
              div(style = "font-size:0.85em;",
                  DT::dataTableOutput("preview_pheno_tbl"))
            )
          ),
          column(6,
            card(
              card_header(tags$b("Pathway Annotation"),
                          tags$small(" (optional; enables enrichment)", style = "color:gray;")),
              p(tags$small("Columns: feature, pathway, annotation")),
              fileInput("file_pathway", NULL, placeholder = "Upload Pathway (.txt/.csv)",
                        accept = c(".txt",".csv")),
              div(class = "d-flex justify-content-between align-items-center mb-1",
                  actionButton("load_example_pathway", "Use example (Human pathways)",
                               icon = icon("database"),
                               class = "btn btn-sm btn-outline-success"),
                  actionButton("reset_pathway", "Reset", icon = icon("times"),
                               class = "btn btn-sm btn-outline-secondary")),
              uiOutput("status_pathway"),
              tags$hr(),
              strong("Preview:"),
              div(style = "font-size:0.85em;",
                  DT::dataTableOutput("preview_pathway_tbl"))
            )
          )
        )
      )
    )
  ),

  # ──────────────────────────────────────────────────────────────────────────
  # Tab 3: Network Files
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "Network Files",
    icon = icon("sitemap"),
    div(
      class = "container-fluid py-3",
      fluidRow(
        column(6,
          card(
            card_header(tags$b("Within-Network"),
                        tags$small(" (edges within Data X)", style = "color:gray;")),
            p(tags$small("Columns: nodeA, nodeB, sign (1 or –1)")),
            fileInput("file_within_net", NULL,
                      placeholder = "Upload within-network (.txt/.csv)",
                      accept = c(".txt",".csv")),
            div(class = "d-flex justify-content-end mb-1",
                actionButton("reset_within_net", "Reset", icon = icon("times"),
                             class = "btn btn-sm btn-outline-secondary")),
            uiOutput("status_within_net"),
            tags$hr(),
            strong("Preview:"),
            div(style = "overflow-x:auto; font-size:0.85em;",
                tableOutput("preview_within_tbl"))
          )
        ),
        column(6,
          uiOutput("between_net_panel")
        )
      ),
      uiOutput("inferred_net_load_ui")
    )
  ),

  # Tab 4: Network Inference
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "Network Inference",
    icon = icon("project-diagram"),
    div(
      class = "container-fluid py-3",
      fluidRow(
        column(4,
          card(
            card_header(tags$b("NetDeconvolute Settings")),
            radioButtons("net_option", "Method:",
                         choices = c("Supervised (glasso only)" = "1",
                                     "Hybrid (glasso + prior)" = "2"),
                         selected = "1"),
            conditionalPanel(
              condition = "input.net_option == '2'",
              fileInput("file_prior_net", "Prior network file:",
                        placeholder = "Upload prior network",
                        accept = c(".txt",".csv")),
              p(tags$small("Required columns: NodeA, NodeB, Dir, NodeA_DT, NodeB_DT"))
            ),
            hr(),
            sliderInput("numLambda", "Number of lambda values:",
                        min = 10, max = 100, value = 30, step = 5),
            checkboxInput("use_lambda_range",
                          "Specify lambda range (optional)",
                          value = FALSE),
            conditionalPanel(
              condition = "input.use_lambda_range",
              fluidRow(
                column(6, numericInput("lambda_min", "λ min:",
                                       value = 0.01, min = 1e-5, step = 0.005)),
                column(6, numericInput("lambda_max", "λ max:",
                                       value = 0.5,  min = 1e-5, step = 0.05))
              ),
              tags$small(class = "text-muted",
                         "Overrides the data-driven lambda grid.")
            ),
            selectInput("criterion", "Model selection criterion:",
                        choices  = c("eBIC","BIC","AIC","CV"),
                        selected = "eBIC"),
            checkboxInput("calibration",
                          "Run calibration (auto-selects lambda)",
                          value = TRUE),
            conditionalPanel(
              condition = "!input.calibration",
              numericInput("optLambda", "Optimal lambda:",
                           value = 0.1, min = 0.001, step = 0.01)
            ),
            hr(),
            checkboxInput("net_logtransform",  "Log-transform input data", value = TRUE),
            checkboxInput("net_standardize",   "Standardize input data",   value = TRUE),
            checkboxInput("net_detectOutliers","Detect outliers (PCA)",    value = TRUE),
            hr(),
            helpText(tags$small(
              icon("info-circle"),
              "Each -omics dataset is z-score standardized per feature independently, then all datasets are row-concatenated into one matrix. The cross-covariance is computed on this combined matrix and graphical lasso is applied to estimate a joint sparse precision matrix."
            )),
            hr(),
            actionButton("run_netdeconv", "Run Network Inference",
                         class = "btn btn-primary w-100",
                         icon  = icon("play"))
          )
        ),
        column(8,
          card(
            card_header(tags$b("Network Inference Log & Results")),
            verbatimTextOutput("net_log"),
            uiOutput("net_calib_ui"),
            uiOutput("net_result_summary"),
            uiOutput("net_heatmaps_ui"),
            uiOutput("net_use_btn"),
            uiOutput("net_download_ui"),
            uiOutput("net_plot_ui")
          )
        )
      )
    )
  ),

  # ──────────────────────────────────────────────────────────────────────────
  # Tab 5: iOmicsPASS Analysis
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "iOmicsPASS Analysis",
    icon = icon("cogs"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        open  = TRUE,
        h5("Analysis Parameters"),
        hr(),
        checkboxInput("do_cv", "Cross-validation", value = TRUE),
        conditionalPanel(
          condition = "input.do_cv",
          sliderInput("num_folds", "Number of CV folds:",
                      min = 2, max = 10, value = 5, step = 1)
        ),
        conditionalPanel(
          condition = "!input.do_cv",
          numericInput("manual_thres", "Min. threshold (skips CV):",
                       value = 2, min = 0, step = 0.1)
        ),
        hr(),
        uiOutput("inferred_net_option_ui"),
        uiOutput("net_enrich_sidebar_ui"),
        tags$hr(style = "margin:6px 0;"),
        uiOutput("gene_enrich_sidebar_ui"),
        hr(),
        checkboxInput("do_knn", "KNN imputation", value = FALSE),
        conditionalPanel(
          condition = "input.do_knn",
          sliderInput("knn_k", "K (neighbours):",
                      min = 2, max = 20, value = 10, step = 1)
        ),
        hr(),
        numericInput("min_obs",  "Min. observations per feature:", value = 5,   min = 1, step = 1),
        numericInput("min_prop", "Min. non-missing proportion:",   value = 0.8, min = 0, max = 1, step = 0.1),
        hr(),
        div(
          class = "d-grid",
          actionButton("run_analysis", "Run Analysis",
                       class = "btn btn-primary btn-lg",
                       icon  = icon("play"))
        ),
        hr(),
        strong("Run log:"),
        verbatimTextOutput("analysis_log")
      ),

      div(
        class = "container-fluid py-3",
        uiOutput("analysis_results_ui")
      )
    )
  ),

  # ──────────────────────────────────────────────────────────────────────────
  # Tab 6: Prediction
  # ──────────────────────────────────────────────────────────────────────────
  nav_panel(
    "Prediction",
    icon = icon("flask"),
    div(
      class = "container py-4",
      style = "max-width:960px;",
      h4("Predict New Samples"),
      p("Upload new omics data for classification using the model produced in the ",
        tags$b("iOmicsPASS Analysis"), " tab."),
      uiOutput("predict_model_status"),
      hr(),
      fluidRow(
        column(4,
          card(
            card_header(tags$b("New Sample Data")),
            p(tags$small(
              "Matrix with features as rows and new samples as columns. ",
              "First column = feature name. ",
              "Must contain all (or most) features from the selected edges."
            )),
            fileInput("file_newdata", "Upload new data:",
                      placeholder = "Upload new data (.txt/.csv)",
                      accept      = c(".txt",".csv")),
            checkboxInput("predict_standardize", "Standardize new data", value = TRUE),
            numericInput("predict_prop", "Min. feature overlap proportion:",
                         value = 0.8, min = 0, max = 1, step = 0.05),
            hr(),
            div(
              class = "d-grid",
              actionButton("run_predict", "Run Prediction",
                           class = "btn btn-primary",
                           icon  = icon("play"))
            )
          )
        ),
        column(8,
          card(
            card_header(tags$b("Prediction Results")),
            uiOutput("predict_results_ui"),
            div(style = "overflow-x:auto;", DTOutput("tbl_predict"))
          )
        )
      )
    )
  ),

)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  rv <- reactiveValues(
    data_X        = NULL,
    data_Y        = NULL,
    data_Z        = NULL,
    phenotype     = NULL,
    pathway       = NULL,
    within_net    = NULL,
    btw_net       = NULL,
    inferred_net  = NULL,
    calib_scores  = NULL,
    cov_mat       = NULL,
    partial_cor_mat = NULL,
    result        = NULL,
    pred_result   = NULL,
    net_log       = "Network inference has not been run yet.\n",
    analysis_log  = "Analysis has not been run yet.\n",
    example_loaded = FALSE,
    edge_classes   = character(0),  # class names from _absdiknew cols
    active_palette      = "classic",  # node palette applied to current network
    active_dir_pal      = "rb",       # edge direction palette applied to current network
    palette_rev         = FALSE,      # whether X/Y node colors are swapped
    net_enrich_result   = NULL,       # R-computed network enrichment (from Run Network Enrichment button)
    gene_enrich_counter = 0L,         # incremented by Run Pathway Enrichment button
    bg_custom           = NULL,       # custom background gene list uploaded by user
    net_export_path     = NULL,       # temp PNG written each time user saves network
    net_layout          = NULL,       # cached igraph layout matrix
    net_layout_key      = NULL,       # node-set hash used to detect layout invalidation
    net_layout_nodes    = NULL        # node order matching net_layout rows
  )

  # ── Load Example Data ──────────────────────────────────────────────────────
  observeEvent(input$load_example, {
    withProgress(message = "Loading example datasets...", value = 0.5, {
      tryCatch({
        e <- new.env(parent = emptyenv())
        data("Tulip_Protein",     package = "iOmicsPASSplus", envir = e)
        data("Tulip_microRNA",    package = "iOmicsPASSplus", envir = e)
        data("PhenotypeFile",     package = "iOmicsPASSplus", envir = e)
        data("PPI_network",       package = "iOmicsPASSplus", envir = e)
        data("TargetScan_network",package = "iOmicsPASSplus", envir = e)
        data("bioPathways",       package = "iOmicsPASSplus", envir = e)

        rv$data_X    <- e$Tulip_Protein
        rv$data_Y    <- e$Tulip_microRNA
        rv$phenotype <- e$PhenotypeFile
        rv$pathway   <- e$bioPathways

        # Standardise network column names to nodeA / nodeB / sign
        ppi <- e$PPI_network
        colnames(ppi) <- c("nodeA", "nodeB", "sign")
        rv$within_net <- ppi

        ts <- e$TargetScan_network
        colnames(ts) <- c("nodeA", "nodeB", "sign")
        rv$btw_net <- ts

        rv$example_loaded <- TRUE

        updateCheckboxInput(session, "use_Y",         value = TRUE)
        updateCheckboxInput(session, "log_X",          value = TRUE)
        updateCheckboxInput(session, "log_Y",          value = TRUE)

        showNotification("Example data loaded. Navigate to each tab to review inputs.",
                         type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Error loading example data:", conditionMessage(e)),
                         type = "error")
      })
    })
  })

  observeEvent(input$load_example_pathway, {
    tryCatch({
      e <- new.env(parent = emptyenv())
      data("bioPathways", package = "iOmicsPASSplus", envir = e)
      rv$pathway <- e$bioPathways
      showNotification("Example Human pathway annotations loaded.", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error loading example pathway:", conditionMessage(e)), type = "error")
    })
  })

  output$example_loaded_msg <- renderUI({
    if (!rv$example_loaded) return(NULL)
    div(class = "alert alert-success mt-3",
        icon("check-circle"), " Example data loaded successfully. ",
        "Navigate to ", tags$b("Network Setup"), " and then ",
        tags$b("iOmicsPASS Analysis"), " to run the workflow.")
  })

  # ── File upload observers ──────────────────────────────────────────────────
  observeEvent(input$file_X, {
    df <- read_upload(input$file_X$datapath)
    if (!is.null(df)) rv$data_X <- df
    else showNotification("Could not parse Data X file.", type = "warning")
  })
  observeEvent(input$file_Y, {
    df <- read_upload(input$file_Y$datapath)
    if (!is.null(df)) rv$data_Y <- df
    else showNotification("Could not parse Data Y file.", type = "warning")
  })
  observeEvent(input$file_Z, {
    df <- read_upload(input$file_Z$datapath)
    if (!is.null(df)) rv$data_Z <- df
    else showNotification("Could not parse Data Z file.", type = "warning")
  })
  observeEvent(input$file_pheno, {
    df <- read_upload(input$file_pheno$datapath)
    if (!is.null(df)) rv$phenotype <- df
    else showNotification("Could not parse Phenotype file.", type = "warning")
  })
  observeEvent(input$file_pathway, {
    df <- read_upload(input$file_pathway$datapath)
    if (is.null(df)) {
      showNotification("Could not parse Pathway file.", type = "warning")
    } else if (ncol(df) < 2L) {
      showNotification("Pathway file must have at least 2 columns (gene, pathwayID).", type = "error")
    } else {
      if (ncol(df) == 2L) df[[3L]] <- df[[2L]]   # use pathwayID as description if column 3 absent
      rv$pathway <- df
    }
  })
  observeEvent(input$file_within_net, {
    df <- read_upload(input$file_within_net$datapath)
    if (is.null(df)) {
      showNotification("Could not parse within-network file.", type = "warning")
    } else if (ncol(df) < 2L) {
      showNotification("Within-network file must have at least 2 columns (nodeA, nodeB).", type = "error")
    } else {
      if (ncol(df) == 2L) {
        # Add default sign column (1 = positive/activating); user can re-upload with signs
        df$sign <- 1L
        showNotification(
          "Within-network has no sign column — defaulting all edges to +1. Re-upload with a 3rd column (1 or -1) if directions are known.",
          type = "warning", duration = 8L
        )
      }
      df[[3L]] <- as.integer(df[[3L]])
      rv$within_net <- df
    }
  })
  observeEvent(input$file_btw_net, {
    df <- read_upload(input$file_btw_net$datapath)
    if (is.null(df)) {
      showNotification("Could not parse between-network file.", type = "warning")
    } else if (ncol(df) < 2L) {
      showNotification("Between-network file must have at least 2 columns (nodeA, nodeB).", type = "error")
    } else {
      if (ncol(df) == 2L) {
        df$sign <- 1L
        showNotification(
          "Between-network has no sign column — defaulting all edges to +1. Re-upload with a 3rd column (1 or -1) if directions are known.",
          type = "warning", duration = 8L
        )
      }
      df[[3L]] <- as.integer(df[[3L]])
      rv$btw_net <- df
    }
  })

  # ── File reset observers ───────────────────────────────────────────────────
  observeEvent(input$reset_X, {
    rv$data_X <- NULL
    session$sendCustomMessage("resetFileInput", "file_X")
  })
  observeEvent(input$reset_Y, {
    rv$data_Y <- NULL
    session$sendCustomMessage("resetFileInput", "file_Y")
  })
  observeEvent(input$reset_Z, {
    rv$data_Z <- NULL
    session$sendCustomMessage("resetFileInput", "file_Z")
  })
  observeEvent(input$reset_pheno, {
    rv$phenotype <- NULL
    session$sendCustomMessage("resetFileInput", "file_pheno")
  })
  observeEvent(input$reset_pathway, {
    rv$pathway <- NULL
    session$sendCustomMessage("resetFileInput", "file_pathway")
  })
  observeEvent(input$reset_within_net, {
    rv$within_net <- NULL
    session$sendCustomMessage("resetFileInput", "file_within_net")
  })
  observeEvent(input$reset_btw_net, {
    rv$btw_net <- NULL
    session$sendCustomMessage("resetFileInput", "file_btw_net")
  })

  # ── Manual lambda ↔ slider sync ───────────────────────────────────────────
  observeEvent(input$manual_lambda_val, {
    req(input$manual_lambda_val, input$netdeconv_lambda)
    if (!isTRUE(all.equal(input$manual_lambda_val, input$netdeconv_lambda, tolerance = 1e-8)))
      updateSliderInput(session, "netdeconv_lambda", value = input$manual_lambda_val)
  }, ignoreInit = TRUE)

  observeEvent(input$netdeconv_lambda, {
    req(input$netdeconv_lambda)
    if (is.null(input$manual_lambda_val) ||
        !isTRUE(all.equal(input$netdeconv_lambda, input$manual_lambda_val, tolerance = 1e-8)))
      updateNumericInput(session, "manual_lambda_val", value = input$netdeconv_lambda)
  }, ignoreInit = TRUE)

  # ── Status indicators ──────────────────────────────────────────────────────
  output$status_X <- renderUI({
    if (!is.null(rv$data_X))
      ok_badge(sprintf("Loaded: %d features x %d samples",
                       nrow(rv$data_X), ncol(rv$data_X) - 1L))
    else
      warn_badge("Required — not loaded")
  })
  output$status_Y <- renderUI({
    if (!isTRUE(input$use_Y)) return(NULL)
    if (!is.null(rv$data_Y))
      ok_badge(sprintf("Loaded: %d features x %d samples",
                       nrow(rv$data_Y), ncol(rv$data_Y) - 1L))
    else
      info_badge("Not loaded")
  })
  output$status_Z <- renderUI({
    if (!isTRUE(input$use_Z)) return(NULL)
    if (!is.null(rv$data_Z))
      ok_badge(sprintf("Loaded: %d features x %d samples",
                       nrow(rv$data_Z), ncol(rv$data_Z) - 1L))
    else
      info_badge("Not loaded")
  })
  output$status_pheno <- renderUI({
    if (!is.null(rv$phenotype)) {
      grp <- tryCatch(table(rv$phenotype[, 2L]), error = function(e) NULL)
      grp_str <- if (!is.null(grp))
        paste(names(grp), grp, sep = "=", collapse = ", ") else ""
      tagList(
        ok_badge(sprintf("Loaded: %d samples", nrow(rv$phenotype))),
        tags$br(),
        tags$small(style = "color:#555;", paste("Groups:", grp_str))
      )
    } else {
      warn_badge("Required — not loaded")
    }
  })
  output$status_pathway <- renderUI({
    if (!is.null(rv$pathway))
      ok_badge(sprintf("Loaded: %d entries", nrow(rv$pathway)))
    else
      info_badge("Not loaded (optional; enables enrichment)")
  })
  output$status_within_net <- renderUI({
    if (!is.null(rv$within_net))
      ok_badge(sprintf("Loaded: %d edges", nrow(rv$within_net)))
    else
      info_badge("Not loaded")
  })
  output$status_btw_net <- renderUI({
    if (!isTRUE(input$use_Y)) return(NULL)
    if (!is.null(rv$btw_net))
      ok_badge(sprintf("Loaded: %d edges", nrow(rv$btw_net)))
    else
      warn_badge("Required when using Data Y")
  })

  # ── Preview tables ─────────────────────────────────────────────────────────
  .preview_dt <- function(df) {
    DT::datatable(df,
                  options  = list(scrollX = TRUE, dom = "t", ordering = FALSE,
                                  pageLength = 5L, columnDefs = list(list(className = "dt-left", targets = "_all"))),
                  rownames = FALSE, class = "table-sm table-striped table-hover")
  }

  output$preview_X_tbl <- DT::renderDataTable({
    req(rv$data_X); .preview_dt(head(rv$data_X, 5L))
  })
  output$preview_Y_tbl <- DT::renderDataTable({
    req(input$use_Y, rv$data_Y); .preview_dt(head(rv$data_Y, 5L))
  })
  output$preview_Z_tbl <- DT::renderDataTable({
    req(input$use_Z, rv$data_Z); .preview_dt(head(rv$data_Z, 5L))
  })
  output$preview_pheno_tbl <- DT::renderDataTable({
    req(rv$phenotype); .preview_dt(head(rv$phenotype, 5L))
  })
  output$preview_pathway_tbl <- DT::renderDataTable({
    req(rv$pathway); .preview_dt(head(rv$pathway, 5L))
  })

  output$preview_within_tbl <- renderTable({
    req(rv$within_net); head(rv$within_net, 5L)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  output$preview_btw_tbl <- renderTable({
    req(rv$btw_net); head(rv$btw_net, 5L)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, rownames = FALSE)

  # ── Between-network panel (conditional on Data Y) ─────────────────────────
  output$between_net_panel <- renderUI({
    if (!isTRUE(input$use_Y)) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 " Enable Data Y in the Data Upload tab to configure cross-omics edges."))
    }
    card(
      card_header(tags$b("Between-Network"),
                  tags$small(" (Data X – Data Y edges; required with Data Y)",
                             style = "color:gray;")),
      p(tags$small("Columns: nodeA (X feature), nodeB (Y feature), sign (1 or –1)")),
      fileInput("file_btw_net", NULL,
                placeholder = "Upload between-network (.txt/.csv)",
                accept = c(".txt",".csv")),
      div(class = "d-flex justify-content-end mb-1",
          actionButton("reset_btw_net", "Reset", icon = icon("times"),
                       class = "btn btn-sm btn-outline-secondary")),
      uiOutput("status_btw_net"),
      tags$hr(),
      strong("Preview:"),
      div(style = "overflow-x:auto; font-size:0.85em;",
          tableOutput("preview_btw_tbl"))
    )
  })

  # ── Network Inference ──────────────────────────────────────────────────────
  observeEvent(input$run_netdeconv, {
    req(rv$data_X)

    rv$net_log      <- "Starting network inference...\n"
    rv$inferred_net <- NULL

    withProgress(message = "Running NetDeconvolute...", value = 0L, {

      tryCatch({
        # Build input list: rows = features, cols = samples
        mat_X <- as.matrix(rv$data_X[, -1L])
        row.names(mat_X) <- rv$data_X[[1L]]
        inputDat <- list(X = mat_X)

        if (isTRUE(input$use_Y) && !is.null(rv$data_Y)) {
          mat_Y <- as.matrix(rv$data_Y[, -1L])
          row.names(mat_Y) <- rv$data_Y[[1L]]
          inputDat[["Y"]] <- mat_Y
        }

        prior_net <- NULL
        if (input$net_option == "2" && !is.null(input$file_prior_net)) {
          prior_net <- read_upload(input$file_prior_net$datapath)
        }

        opt_lam <- if (isTRUE(input$calibration)) NULL else input$optLambda

        lam_vec <- if (isTRUE(input$use_lambda_range) &&
                       !is.na(input$lambda_min) && !is.na(input$lambda_max) &&
                       input$lambda_min < input$lambda_max)
          exp(seq(log(input$lambda_max), log(input$lambda_min),
                  length.out = as.integer(input$numLambda)))
        else
          NULL

        log_msgs <- character(0L)

        old_wd <- getwd()
        setwd(app_dir)
        on.exit(setwd(old_wd), add = TRUE)

        withCallingHandlers({
          incProgress(0.2, detail = "Estimating network...")
          NetDeconvolute(
            inputDat        = inputDat,
            option          = as.integer(input$net_option),
            NetworkFile     = prior_net,
            Calibration     = isTRUE(input$calibration),
            optLambda       = opt_lam,
            lambda.vec      = lam_vec,
            criterion       = input$criterion,
            log.transform   = isTRUE(input$net_logtransform),
            standardization = isTRUE(input$net_standardize),
            detectOutliers  = isTRUE(input$net_detectOutliers),
            numLambda       = as.integer(input$numLambda),
            Plot.PCA        = TRUE,
            Plot.modelSelect = TRUE,
            Plot.Covariance = FALSE,
            verbose         = TRUE
          )
          incProgress(0.9, detail = "Done.")
        }, message = function(m) {
          log_msgs <<- c(log_msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        })

        rv$net_log <- paste(c(log_msgs,
                              "\n--- Network inference completed. ---\n"),
                            collapse = "")

        # Try to read inferred network file
        net_file <- if (input$net_option == "2")
          file.path(app_dir, "Output", "Fusednetwork_hybridmethod.txt")
        else
          file.path(app_dir, "Output", "Estimated_Network_glasso.txt")

        if (file.exists(net_file)) {
          rv$inferred_net <- tryCatch(
            read.delim(net_file, as.is = TRUE),
            error = function(e) NULL
          )
        }

        calib_file <- file.path(app_dir, "Output", "CalibrationScores_glasso.csv")
        if (file.exists(calib_file)) {
          rv$calib_scores <- tryCatch(
            read.csv(calib_file, as.is = TRUE),
            error = function(e) NULL
          )
        }

        showNotification("Network inference complete!", type = "message", duration = 3L)

      }, error = function(e) {
        rv$net_log <- paste0(rv$net_log, "\nERROR: ", conditionMessage(e), "\n")
        showNotification(paste("Network inference error:", conditionMessage(e)),
                         type = "error")
      })
    })
  })

  output$net_log <- renderText({ rv$net_log })

  output$net_result_summary <- renderUI({
    if (is.null(rv$inferred_net)) return(NULL)
    net      <- rv$inferred_net
    nodes_a  <- as.character(net[[1L]])
    nodes_b  <- as.character(net[[2L]])
    n_e      <- nrow(net)
    all_nodes <- unique(c(nodes_a, nodes_b))
    n_n      <- length(all_nodes)
    max_e    <- n_n * (n_n - 1L) / 2L
    sparsity <- round(100 * (1 - n_e / max_e), 1)

    # EdgeType column is written by NetDeconvolute; fall back to node-name heuristic
    if ("EdgeType" %in% colnames(net)) {
      n_within  <- sum(net$EdgeType == "Within",  na.rm = TRUE)
      n_between <- sum(net$EdgeType == "Between", na.rm = TRUE)
    } else {
      is_x <- function(v) grepl("_dataX$", v, ignore.case = TRUE)
      is_y <- function(v) grepl("_dataY$", v, ignore.case = TRUE)
      n_within  <- sum((is_x(nodes_a) & is_x(nodes_b)) |
                       (is_y(nodes_a) & is_y(nodes_b)))
      n_between <- n_e - n_within
    }
    pct_w <- round(100 * n_within  / max(n_e, 1L), 1)
    pct_b <- round(100 * n_between / max(n_e, 1L), 1)

    tagList(
      div(class = "alert alert-success mt-2",
          icon("check-circle"),
          sprintf(" Inferred network: %d edges, %d unique nodes.", n_e, n_n)),
      div(class = "table-responsive mt-1 mb-2",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$thead(tags$tr(
              tags$th("Metric"), tags$th("Value")
            )),
            tags$tbody(
              tags$tr(tags$td("Total edges"),    tags$td(n_e)),
              tags$tr(tags$td("Within-network edges"),
                      tags$td(sprintf("%d  (%.1f%%)", n_within, pct_w))),
              tags$tr(tags$td("Between-network edges"),
                      tags$td(sprintf("%d  (%.1f%%)", n_between, pct_b))),
              tags$tr(tags$td("Max possible edges (p(p−1)/2)"), tags$td(max_e)),
              tags$tr(tags$td("Sparsity"),
                      tags$td(sprintf("%.1f%%", sparsity)))
            )
          )),
      tags$b("First 10 rows of estimated network:"),
      DT::dataTableOutput("net_preview_table")
    )
  })

  output$net_preview_table <- DT::renderDataTable({
    req(rv$inferred_net)
    DT::datatable(head(rv$inferred_net, 10L),
                  options = list(dom = "t", scrollX = TRUE, ordering = FALSE),
                  rownames = FALSE,
                  class = "table-sm table-striped")
  })

  output$net_use_btn <- renderUI({
    if (is.null(rv$inferred_net)) return(NULL)
    div(class = "mt-2",
        tags$small(class = "text-muted",
                   icon("info-circle"),
                   " Clicking below will: set the inferred network as within-network, load the combined data as Data X, and clear Data Y and the between-network file."),
        tags$br(),
        actionButton("use_inferred_net",
                     "Use inferred network + load combined data as Data X →",
                     class = "btn btn-success mt-1 w-100",
                     icon  = icon("arrow-right")))
  })

  output$net_download_ui <- renderUI({
    net_file <- file.path(app_dir, "Output", "Estimated_Network_glasso.txt")
    if (!file.exists(net_file)) return(NULL)
    div(class = "mt-2",
        downloadButton("download_net", "Download Estimated Network",
                       class = "btn btn-outline-secondary btn-sm",
                       icon  = icon("download")))
  })

  output$download_net <- downloadHandler(
    filename = "Estimated_Network_glasso.txt",
    content  = function(file) {
      net_file <- file.path(app_dir, "Output", "Estimated_Network_glasso.txt")
      file.copy(net_file, file)
    }
  )

  output$inferred_net_load_ui <- renderUI({
    if (is.null(rv$inferred_net)) return(NULL)
    net <- rv$inferred_net
    n_edges <- nrow(net)
    already_loaded <- !is.null(rv$within_net) &&
      identical(colnames(rv$within_net)[1L:2L], c("nodeA", "nodeB"))
    tagList(
      tags$hr(),
      card(
        card_header(
          tags$b(icon("check-circle", style = "color:#27ae60;"),
                 " Inferred Network Available")
        ),
        p(tags$small(
          sprintf("Network with %d edges estimated in the Network Inference tab.", n_edges),
          tags$br(),
          "Clicking the button below will set the inferred network as within-network,",
          " load the combined (X+Y) data as Data X, and clear Data Y and the between-network file."
        )),
        if (already_loaded)
          p(tags$small(class = "text-success",
                       icon("check"), " Already applied as within-network."))
        else
          NULL,
        actionButton("load_inferred_to_files",
                     if (already_loaded) "Re-apply inferred network + combined data"
                     else "Apply inferred network + load combined data as Data X",
                     class = "btn btn-outline-success btn-sm",
                     icon  = icon("upload"))
      )
    )
  })

  # Shared helper: push inferred network to within_net + load combined data as data_X
  .apply_inferred_net <- function() {
    net <- rv$inferred_net
    if (is.null(net)) return()

    # Set within-network (cols 1-3: nodeA, nodeB, dir/sign)
    if (ncol(net) >= 3L) {
      tmp <- net[, 1L:3L]
      colnames(tmp) <- c("nodeA", "nodeB", "sign")
      tmp[[3L]] <- as.integer(tmp[[3L]])
      # Remove edges with direction 0 (sign(pp) = 0 for near-zero partial correlations)
      tmp <- tmp[tmp[[3L]] != 0L, , drop = FALSE]
      rv$within_net <- tmp
    }

    # Load Combined_data.txt (X+Y concatenated & standardised) as the new Data X
    combined_file <- file.path(app_dir, "Output", "Combined_data.txt")
    if (file.exists(combined_file)) {
      combined_df <- tryCatch(
        read.delim(combined_file, as.is = TRUE, check.names = FALSE),
        error = function(e) NULL
      )
      if (!is.null(combined_df)) {
        rv$data_X  <- combined_df
        rv$data_Y  <- NULL
        rv$btw_net <- NULL
        updateCheckboxInput(session, "use_Y", value = FALSE)
        showNotification(
          sprintf(
            "Combined data (%d features, %d samples) loaded as Data X. Data Y and between-network cleared. Proceed to iOmicsPASS Analysis.",
            nrow(combined_df), ncol(combined_df) - 1L
          ),
          type = "message", duration = 6L
        )
        return()
      }
    }
    showNotification(
      "Inferred network set as within-network. (Combined_data.txt not found — Data X unchanged.)",
      type = "warning", duration = 5L
    )
  }

  observeEvent(input$load_inferred_to_files, {
    req(rv$inferred_net)
    .apply_inferred_net()
  })

  output$net_plot_ui <- renderUI({
    png_file <- file.path(app_dir, "Output", "Plots_glasso.png")
    if (!file.exists(png_file)) return(NULL)
    tagList(
      tags$hr(),
      tags$b("glasso Network Plot:"),
      tags$br(),
      tags$img(src = png_file, style = "max-width:100%; margin-top:8px;")
    )
  })

  observeEvent(input$use_inferred_net, {
    req(rv$inferred_net)
    .apply_inferred_net()
  })

  observeEvent(input$use_inferred_analysis, {
    req(rv$inferred_net)
    .apply_inferred_net()
  })

  # ── Calibration criteria plots + lambda selector ──────────────────────────
  output$net_calib_ui <- renderUI({
    req(rv$calib_scores)
    df   <- rv$calib_scores
    lams <- df$lambda
    lam_aic  <- lams[which.min(df$AIC)]
    lam_bic  <- lams[which.min(df$BIC)]
    lam_ebic <- lams[which.min(df$eBIC)]
    lam_cv   <- lams[which.max(df$CV)]
    def_lam  <- switch(isolate(input$criterion),
                       AIC  = lam_aic,  BIC = lam_bic,
                       eBIC = lam_ebic, CV  = lam_cv,
                       lam_ebic)
    cur_lam <- if (!is.null(isolate(input$netdeconv_lambda))) isolate(input$netdeconv_lambda) else def_lam
    tagList(
      hr(),
      tags$b("Model Selection Criteria — select an optimal λ"),
      plotOutput("net_calib_plot", height = "400px"),
      fluidRow(
        column(6,
          sliderInput("netdeconv_lambda", "Slide to select λ:",
                      min   = round(min(lams, na.rm = TRUE), 5),
                      max   = round(max(lams, na.rm = TRUE), 5),
                      value = round(cur_lam, 5),
                      step  = round((max(lams, na.rm = TRUE) - min(lams, na.rm = TRUE)) / 200, 5),
                      width = "100%")
        ),
        column(3,
          numericInput("manual_lambda_val", "Or enter λ manually:",
                       value = round(cur_lam, 5),
                       min   = round(min(lams, na.rm = TRUE), 5),
                       max   = round(max(lams, na.rm = TRUE), 5),
                       step  = round((max(lams, na.rm = TRUE) - min(lams, na.rm = TRUE)) / 200, 5))
        ),
        column(3,
          div(class = "mt-2",
              tags$small(tags$b("Suggested λ:"), tags$br(),
                         sprintf("AIC:  %.4f", lam_aic),  tags$br(),
                         sprintf("BIC:  %.4f", lam_bic),  tags$br(),
                         sprintf("eBIC: %.4f", lam_ebic), tags$br(),
                         sprintf("CV:   %.4f", lam_cv))),
          div(class = "d-grid mt-2",
              actionButton("run_with_lambda",
                           "Run inference with selected λ",
                           class = "btn btn-primary btn-sm",
                           icon  = icon("play")))
        )
      )
    )
  })

  output$net_calib_plot <- renderPlot({
    req(rv$calib_scores)
    df  <- rv$calib_scores
    lam <- if (!is.null(input$netdeconv_lambda)) input$netdeconv_lambda else NA
    par(mfrow = c(2L, 2L), mar = c(3, 3, 2.5, 1), mgp = c(1.8, 0.5, 0))
    .draw_crit <- function(x, y, main, best_fn) {
      ok <- !is.na(y)
      plot(x[ok], y[ok], type = "o", pch = 16, cex = 0.7, col = "#7f8c8d",
           xlab = "", ylab = "", main = main, axes = TRUE)
      bi <- best_fn(y)
      points(x[bi], y[bi], pch = 16, cex = 1.6, col = "#3498db")
      if (!is.na(lam)) abline(v = lam, col = "#e74c3c", lty = 2, lwd = 1.5)
    }
    .draw_crit(df$lambda, df$AIC,  "AIC",          which.min)
    .draw_crit(df$lambda, df$BIC,  "BIC",          which.min)
    .draw_crit(df$lambda, df$eBIC, "extended-BIC", which.min)
    .draw_crit(df$lambda, df$CV,   "k-fold CV",    which.max)
  })

  observeEvent(input$run_with_lambda, {
    req(rv$data_X, rv$calib_scores, input$netdeconv_lambda)
    lam <- input$netdeconv_lambda
    rv$net_log <- paste0(rv$net_log,
                         sprintf("\nRunning inference with λ = %.5f...\n", lam))
    rv$inferred_net  <- NULL
    rv$cov_mat       <- NULL
    rv$partial_cor_mat <- NULL

    withProgress(message = "Running NetDeconvolute...", value = 0L, {
      tryCatch({
        mat_X <- as.matrix(rv$data_X[, -1L])
        row.names(mat_X) <- rv$data_X[[1L]]
        inputDat <- list(X = mat_X)
        if (isTRUE(input$use_Y) && !is.null(rv$data_Y)) {
          mat_Y <- as.matrix(rv$data_Y[, -1L])
          row.names(mat_Y) <- rv$data_Y[[1L]]
          inputDat[["Y"]] <- mat_Y
        }
        prior_net <- NULL
        if (input$net_option == "2" && !is.null(input$file_prior_net)) {
          prior_net <- read_upload(input$file_prior_net$datapath)
        }
        log_msgs <- character(0L)
        old_wd <- getwd()
        setwd(app_dir)
        on.exit(setwd(old_wd), add = TRUE)
        withCallingHandlers({
          incProgress(0.5, detail = "Estimating network...")
          NetDeconvolute(
            inputDat        = inputDat,
            option          = as.integer(input$net_option),
            NetworkFile     = prior_net,
            Calibration     = FALSE,
            optLambda       = lam,
            criterion       = input$criterion,
            log.transform   = isTRUE(input$net_logtransform),
            standardization = isTRUE(input$net_standardize),
            detectOutliers  = FALSE,
            numLambda       = as.integer(input$numLambda),
            Plot.PCA        = FALSE,
            Plot.modelSelect = FALSE,
            Plot.Covariance = FALSE,
            verbose         = TRUE
          )
          incProgress(0.9, detail = "Done.")
        }, message = function(m) {
          log_msgs <<- c(log_msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        })
        rv$net_log <- paste(c(rv$net_log, log_msgs,
                               "\n--- Network inference completed. ---\n"),
                             collapse = "")
        net_file <- if (input$net_option == "2")
          file.path(app_dir, "Output", "Fusednetwork_hybridmethod.txt")
        else
          file.path(app_dir, "Output", "Estimated_Network_glasso.txt")
        if (file.exists(net_file)) {
          rv$inferred_net <- tryCatch(
            read.delim(net_file, as.is = TRUE), error = function(e) NULL)
        }

        # Read partial correlation matrix
        pc_file <- file.path(app_dir, "Output", "PartialCorrelation_icov.txt")
        if (file.exists(pc_file)) {
          rv$partial_cor_mat <- tryCatch(
            read.delim(pc_file, check.names = FALSE, as.is = TRUE),
            error = function(e) NULL)
        }

        # Compute cross-covariance matrix from standardised data (mirrors NetDeconvolute)
        tryCatch({
          .std <- function(mat, logt, stdt) {
            if (logt) mat <- log(mat)
            if (stdt) {
              mn <- rowMeans(mat, na.rm = TRUE)
              sd_v <- apply(mat, 1L, sd, na.rm = TRUE)
              mat <- sweep(mat, 1L, mn) / sd_v
            }
            mat
          }
          mx <- .std(as.matrix(rv$data_X[, -1L]), isTRUE(input$net_logtransform),
                     isTRUE(input$net_standardize))
          rownames(mx) <- rv$data_X[[1L]]
          combined <- mx
          if (isTRUE(input$use_Y) && !is.null(rv$data_Y)) {
            my <- .std(as.matrix(rv$data_Y[, -1L]), isTRUE(input$net_logtransform),
                       isTRUE(input$net_standardize))
            rownames(my) <- rv$data_Y[[1L]]
            combined <- rbind(combined, my)
          }
          rv$cov_mat <- cov(t(combined), use = "pairwise.complete.obs")
        }, error = function(e) NULL)

        showNotification("Network inference complete!", type = "message", duration = 3L)
      }, error = function(e) {
        rv$net_log <- paste0(rv$net_log, "\nERROR: ", conditionMessage(e), "\n")
        showNotification(paste("Error:", conditionMessage(e)), type = "error")
      })
    })
  })

  # ── Covariance and precision matrix heatmaps ──────────────────────────────
  output$net_heatmaps_ui <- renderUI({
    if (is.null(rv$cov_mat) && is.null(rv$partial_cor_mat)) return(NULL)
    tagList(
      hr(),
      tags$b("Matrix Heatmaps"),
      fluidRow(
        column(6,
               tags$b("Cross-covariance matrix"), tags$br(),
               plotOutput("net_cov_plot", height = "420px")),
        column(6,
               tags$b("Partial correlation (precision) matrix"), tags$br(),
               plotOutput("net_prec_plot", height = "420px"))
      )
    )
  })

  output$net_cov_plot <- renderPlot({
    req(rv$cov_mat)
    mat <- rv$cov_mat
    if (nrow(mat) > 100L) {
      top <- order(diag(mat), decreasing = TRUE)[1:100]
      mat <- mat[top, top]
    }
    p <- pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                            main = "Cross-covariance",
                            show_rownames = nrow(mat) <= 40L,
                            show_colnames = ncol(mat) <= 40L,
                            silent = TRUE)
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })

  output$net_prec_plot <- renderPlot({
    req(rv$partial_cor_mat)
    pc  <- rv$partial_cor_mat
    rn  <- as.character(pc[[2L]])
    mat <- as.matrix(pc[, 3:ncol(pc)])
    rownames(mat) <- rn
    colnames(mat) <- rn
    diag(mat) <- 0
    if (nrow(mat) > 100L) {
      top <- order(rowSums(abs(mat)), decreasing = TRUE)[1:100]
      mat <- mat[top, top]
    }
    p <- pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                            main = "Partial correlations",
                            show_rownames = nrow(mat) <= 40L,
                            show_colnames = ncol(mat) <= 40L,
                            silent = TRUE)
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })

  # ── Analysis sidebar: use inferred network option ─────────────────────────
  output$inferred_net_option_ui <- renderUI({
    if (is.null(rv$inferred_net)) return(NULL)
    already <- !is.null(rv$within_net)
    div(class = "alert alert-info p-2 mb-1",
        tags$small(icon("project-diagram"), tags$b(" Inferred network available")),
        tags$br(),
        if (!already)
          tagList(
            tags$small("Applies inferred network as within-network and loads combined data as Data X."),
            div(class = "d-grid mt-1",
                actionButton("use_inferred_analysis", "Apply inferred network + combined data",
                             class = "btn btn-sm btn-outline-primary",
                             icon  = icon("check")))
          )
        else
          tags$small(tags$em(class = "text-success", icon("check-circle"),
                             " Applied — combined data loaded as Data X"))
    )
  })

  # ── Enrichment sidebar buttons ─────────────────────────────────────────────
  output$net_enrich_sidebar_ui <- renderUI({
    has_pw     <- !is.null(rv$pathway)
    has_result <- !is.null(rv$result) && !is.null(rv$result$EdgesSelected)
    note <- if (!has_pw) " — load pathway file"
            else if (!has_result) " — run analysis first"
            else NULL
    tagList(
      tags$small(tags$b("Network Enrichment"),
                 if (!is.null(note)) tags$span(note, style = "color:gray;") else NULL),
      div(class = "d-grid mt-1",
          actionButton("run_net_enrich", "Run Network Enrichment",
                       class = if (!has_pw || !has_result)
                         "btn btn-outline-secondary btn-sm disabled"
                       else "btn btn-outline-primary btn-sm",
                       icon = icon("project-diagram")))
    )
  })

  output$gene_enrich_sidebar_ui <- renderUI({
    has_pw     <- !is.null(rv$pathway)
    has_result <- !is.null(rv$result) && !is.null(rv$result$EdgesSelected)
    note <- if (!has_pw) " — load pathway file"
            else if (!has_result) " — run analysis first"
            else NULL
    tagList(
      tags$small(tags$b("Gene-level Enrichment"),
                 if (!is.null(note)) tags$span(note, style = "color:gray;") else NULL),
      tags$small("Background:"),
      radioButtons("gene_enrich_bg", NULL,
                   choices  = c("Data features"                          = "data",
                                "All pathway genes"                      = "pathway",
                                "All HGNC genes (transcriptomics)"       = "hgnc_all",
                                "Genes with protein product (proteomics)" = "hgnc_protein_product",
                                "Custom list (upload)"                   = "custom"),
                   selected = "data", inline = FALSE),
      conditionalPanel(
        condition = "input.gene_enrich_bg === 'custom'",
        fileInput("bg_genes_file", NULL,
                  accept      = c(".txt", ".csv", ".tsv"),
                  placeholder = "One gene per line",
                  width       = "100%"),
        tags$small(class = "text-muted",
                   "Plain text file, one gene name per line (or first column of a table).")
      ),
      div(class = "d-grid mt-1",
          actionButton("run_gene_enrich", "Run Pathway Enrichment",
                       class = if (!has_pw || !has_result)
                         "btn btn-outline-secondary btn-sm disabled"
                       else "btn btn-outline-success btn-sm",
                       icon = icon("dna")))
    )
  })

  observeEvent(input$bg_genes_file, {
    req(input$bg_genes_file)
    tryCatch({
      lines <- readLines(input$bg_genes_file$datapath, warn = FALSE)
      lines <- trimws(lines[nzchar(trimws(lines))])
      # If first line looks like a header with tab, treat as table and take col 1
      if (any(grepl("\t", lines))) {
        tbl   <- read.table(input$bg_genes_file$datapath, sep = "\t",
                            header = TRUE, stringsAsFactors = FALSE, quote = "")
        lines <- unique(as.character(tbl[[1L]]))
      }
      rv$bg_custom <- unique(lines)
      showNotification(sprintf("Custom background: %d genes loaded.", length(rv$bg_custom)),
                       type = "message", duration = 4L)
    }, error = function(e) {
      showNotification(paste("Error reading background file:", conditionMessage(e)),
                       type = "error", duration = 6L)
    })
  })

  # ── Run Network Enrichment (R implementation) ──────────────────────────────
  observeEvent(input$run_net_enrich, {
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) {
      showNotification("Run analysis first.", type = "warning"); return()
    }
    if (is.null(rv$pathway)) {
      showNotification("Load a pathway file first.", type = "warning"); return()
    }

    withProgress(message = "Computing network enrichment...", value = 0, {
      tryCatch({
        es      <- rv$result$EdgesSelected
        pathway <- rv$pathway
        correct <- input$net_enrich_correct %||% "none"

        # Background edges from both networks
        build_edges <- function(net) {
          if (is.null(net) || nrow(net) == 0L) return(character(0L))
          a <- as.character(net[[1L]]); b <- as.character(net[[2L]])
          unique(paste(pmin(a, b), pmax(a, b), sep = "___"))
        }
        all_edges <- unique(c(build_edges(rv$within_net), build_edges(rv$btw_net)))
        N <- length(all_edges)
        if (N == 0L) {
          showNotification("No background edges found in loaded networks.", type = "error")
          return()
        }

        pw_gene <- as.character(pathway[[1L]])
        pw_id   <- as.character(pathway[[2L]])
        pw_name <- as.character(pathway[[3L]])
        pw_list  <- split(pw_gene, pw_id)
        pw_names <- tapply(pw_name, pw_id, function(x) x[[1L]])

        # Edge-oriented pathway map
        edge_a <- sub("___.*", "", all_edges)
        edge_b <- sub(".*___", "", all_edges)
        pw_edge_map   <- lapply(pw_list, function(g) {
          g <- unique(g)
          all_edges[edge_a %in% g & edge_b %in% g]
        })
        pw_edge_sizes <- lengths(pw_edge_map)

        grp_cols <- grep("(?i)_absdiknew$", names(es), value = TRUE, perl = TRUE)
        groups   <- sub("(?i)_absdiknew$", "", grp_cols, perl = TRUE)

        incProgress(0.3, detail = "Running hypergeometric tests...")

        result_up   <- vector("list", length(groups)); names(result_up)   <- groups
        result_down <- vector("list", length(groups)); names(result_down) <- groups
        clean_node  <- function(x) sub("\\.\\.\\.(datX|datY)$", "", x)

        for (grp in groups) {
          dir_col <- paste0(grp, "_dir")
          for (direction in c("up", "down")) {
            sub_es <- if (dir_col %in% names(es))
              es[!is.na(es[[dir_col]]) & es[[dir_col]] == direction, , drop = FALSE]
            else es
            if (nrow(sub_es) == 0L) next

            nodeA <- clean_node(sub_es$NodeA); nodeB <- clean_node(sub_es$NodeB)
            sig_edges <- unique(paste(pmin(nodeA, nodeB), pmax(nodeA, nodeB), sep = "___"))
            sig_edges <- sig_edges[sig_edges %in% all_edges]
            n <- length(sig_edges)
            if (n == 0L) next

            np   <- length(pw_list)
            p_v  <- numeric(np); k_v <- integer(np); gene_str <- character(np)
            valid <- logical(np)

            for (i in seq_len(np)) {
              K <- pw_edge_sizes[[i]]
              if (K < 3L) next
              valid[i] <- TRUE
              k <- sum(sig_edges %in% pw_edge_map[[i]])
              k_v[i] <- k
              p_v[i] <- if (k == 0L) 1.0 else phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
              if (k > 0L) {
                sp <- sig_edges[sig_edges %in% pw_edge_map[[i]]]
                gs <- unique(c(sub("___.*", "", sp), sub(".*___", "", sp)))
                gene_str[i] <- paste(sort(intersect(gs, unique(pw_list[[i]]))), collapse = ";")
              }
            }

            p_adj <- p_v
            vi    <- which(valid)
            if (correct != "none" && length(vi) > 1L) {
              p_adj[vi] <- if (correct == "qvalue") {
                tryCatch(
                  qvalue::qvalue(p_v[vi], pi0.method = "smoother")$qvalues,
                  error = function(e) p.adjust(p_v[vi], method = "BH")
                )
              } else {
                p.adjust(p_v[vi], method = correct)
              }
            }

            df <- data.frame(
              Pathway           = names(pw_list),
              PathAnnotation    = as.character(pw_names[names(pw_list)]),
              Direction         = direction,
              HypergeoPval      = p_v,
              AdjPval           = p_adj,
              Enriched_Edgesize = k_v,
              NumGene_inPathway = pw_edge_sizes,
              Geneset_Enriched  = gene_str,
              stringsAsFactors  = FALSE
            )[valid, , drop = FALSE]
            df <- df[order(df$HypergeoPval), , drop = FALSE]

            if (direction == "up") result_up[[grp]]   <- df
            else                   result_down[[grp]] <- df
          }
        }

        rv$net_enrich_result <- list(Enrichment_up = result_up, Enrichment_down = result_down)
        showNotification("Network enrichment complete.", type = "message")

      }, error = function(e) {
        showNotification(paste("Network enrichment failed:", conditionMessage(e)), type = "error")
      })
    })
  })

  # ── Run Gene-level Enrichment (button trigger) ─────────────────────────────
  observeEvent(input$run_gene_enrich, {
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) {
      showNotification("Run analysis first.", type = "warning"); return()
    }
    if (is.null(rv$pathway)) {
      showNotification("Load a pathway file first.", type = "warning"); return()
    }
    rv$gene_enrich_counter <- rv$gene_enrich_counter + 1L
  })

  # ── Run iOmicsPASS Analysis ────────────────────────────────────────────────
  observeEvent(input$run_analysis, {

    # Input validation
    errs <- character(0L)
    if (is.null(rv$data_X))                               errs <- c(errs, "Data X is required.")
    if (is.null(rv$phenotype))                            errs <- c(errs, "Phenotype file is required.")
    if (is.null(rv$within_net) && is.null(rv$btw_net))   errs <- c(errs, "At least one network (within or between) is required.")
    if (!is.null(rv$within_net) && ncol(rv$within_net) < 3L)
      errs <- c(errs, "Within-network must have 3 columns (nodeA, nodeB, sign). Please re-upload.")
    if (!is.null(rv$btw_net) && ncol(rv$btw_net) < 3L)
      errs <- c(errs, "Between-network must have 3 columns (nodeA, nodeB, sign). Please re-upload.")
    if (is.null(rv$phenotype) || ncol(rv$phenotype) < 2L)
      errs <- c(errs, "Phenotype file must have at least 2 columns (SampleID, Group).")
    if (isTRUE(input$use_Y) && is.null(rv$data_Y))       errs <- c(errs, "Data Y is enabled but not loaded.")
    if (isTRUE(input$use_Y) && is.null(rv$btw_net))      errs <- c(errs, "Between-network is required when using Data Y.")

    if (length(errs) > 0L) {
      showNotification(paste(errs, collapse = "\n"), type = "error", duration = 8L)
      return()
    }

    # Keep old result visible during analysis (avoids double re-render cycle)
    rv$analysis_log      <- "Preparing inputs...\n"

    withProgress(message = "Running iOmicsPASS analysis...", value = 0L, {

      tryCatch({

        data_Y_use    <- if (isTRUE(input$use_Y))  rv$data_Y  else NULL
        data_Z_use    <- if (isTRUE(input$use_Z))  rv$data_Z  else NULL
        btw_net_use   <- if (isTRUE(input$use_Y))  rv$btw_net else NULL
        norm_Z_use    <- if (isTRUE(input$use_Z))  input$normalizeZby else NULL
        # Enrichment is handled separately via the Run Network/Pathway Enrichment buttons
        do_enrich   <- FALSE
        pathway_use <- NULL
        thres_use     <- if (!isTRUE(input$do_cv)) input$manual_thres else NULL
        log_tf        <- c(isTRUE(input$log_X),
                           isTRUE(input$log_Y),
                           isTRUE(input$log_Z))
        knn_k_use     <- if (isTRUE(input$do_knn)) as.integer(input$knn_k) else 10L

        log_msgs <- character(0L)

        old_wd <- getwd()
        setwd(app_dir)
        on.exit(setwd(old_wd), add = TRUE)

        incProgress(0.1, detail = "Creating input parameter file...")

        withCallingHandlers({
          createInputParam(
            data.X          = rv$data_X,
            data.Y          = data_Y_use,
            data.Z          = data_Z_use,
            phenotype       = rv$phenotype,
            within.net      = rv$within_net,
            btw.net         = btw_net_use,
            pathway         = pathway_use,
            standardize.data = isTRUE(input$std_data),
            log.transform   = log_tf,
            normalizeBy     = norm_Z_use,
            Cross.Validate  = isTRUE(input$do_cv),
            num.Kfold       = as.integer(input$num_folds),
            min.thres       = thres_use,
            Enrichment      = do_enrich,
            knn.impute      = isTRUE(input$do_knn),
            knn.k           = knn_k_use,
            min.obs         = as.integer(input$min_obs),
            min.prop        = input$min_prop,
            tag             = "gui",
            return.list     = FALSE
          )
        }, message = function(m) {
          log_msgs <<- c(log_msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        })

        incProgress(0.3, detail = "Running iOmicsPASS core (may take a few minutes)...")

        withCallingHandlers({
          res <- iOmicsPASS.R(
            ff             = "input_param_gui",
            Cross.Validate = isTRUE(input$do_cv),
            plotCV         = FALSE,
            use.rcpp       = TRUE
          )
          rv$result <- res
          rv$net_enrich_result <- NULL   # clear stale enrichment alongside new result
        }, message = function(m) {
          log_msgs <<- c(log_msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        })

        incProgress(1.0, detail = "Complete.")
        rv$analysis_log <- paste(c(log_msgs,
                                   "\n--- Analysis completed successfully. ---\n"),
                                 collapse = "")
        showNotification("iOmicsPASS analysis complete!", type = "message", duration = 3L)

      }, error = function(e) {
        rv$analysis_log <- paste0(rv$analysis_log, "\nERROR: ", conditionMessage(e), "\n")
        showNotification(paste("Analysis error:", conditionMessage(e)),
                         type = "error", duration = 10L)
      })
    })
  })

  output$analysis_log <- renderText({ rv$analysis_log })

  # ── Analysis Results UI container ─────────────────────────────────────────
  output$analysis_results_ui <- renderUI({
    if (is.null(rv$result)) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 " Set parameters on the left and click ",
                 tags$b("Run Analysis"),
                 " to see results here."))
    }

    navset_tab(
      nav_panel(
        tags$span("CV Plot", icon("chart-line", style = "margin-left:4px;")),
        div(class = "mt-3",
            uiOutput("cv_plot_or_msg"))
      ),
      nav_panel(
        "Selected Edges",
        div(class = "mt-2",
            div(class = "d-flex align-items-center gap-2 mb-1",
                downloadButton("dl_edges", "Download edges table", class = "btn btn-sm btn-outline-primary"),
                textOutput("edges_count_badge", inline = TRUE)),
            div(class = "mt-2", DTOutput("tbl_edges")))
      ),
      nav_panel(
        "Sample Classification",
        div(class = "mt-2",
            div(class = "mb-2",
                downloadButton("dl_sample_class", "Download classification table (.txt)",
                               class = "btn btn-sm btn-outline-primary",
                               icon  = icon("download"))),
            uiOutput("confusion_matrix_ui"),
            div(class = "mt-2", DTOutput("tbl_probs")))
      ),
      nav_panel(
        "Network Summary",
        div(class = "mt-3",
            div(class = "d-flex flex-wrap gap-3 mb-2",
                uiOutput("net_group_select_ui"),
                div(
                  selectInput("net_palette", "Node color palette:",
                              choices = c("Classic (Blue/Orange)"  = "classic",
                                          "Blue/Purple"            = "purple",
                                          "Cool (Teal/Coral)"      = "cool",
                                          "Monochrome (Dark/Light)" = "mono"),
                              width = "220px"),
                  div(class = "d-flex gap-2 mt-1",
                      actionButton("apply_palette", "Apply palette",
                                   class = "btn btn-sm btn-outline-secondary",
                                   icon  = icon("palette")),
                      actionButton("reverse_palette", "Swap colors",
                                   class = "btn btn-sm btn-outline-secondary",
                                   icon  = icon("right-left")))),
                selectInput("edge_dir_pal", "Edge direction colors:",
                            choices = c("Red (up) / Blue (down)"       = "rb",
                                        "Blue (up) / Red (down)"       = "br",
                                        "Green (up) / Red (down)"      = "gr",
                                        "Blue (up) / Gold (down)"      = "bg",
                                        "Teal (up) / Blue (down)"      = "tg"),
                            width = "220px"),
                selectInput("edge_smooth_type", "Edge style:",
                            choices = c("Curved" = "continuous", "Straight" = "none"),
                            width = "140px"),
                uiOutput("edge_sign_toggle_ui"),
                div(class = "d-flex gap-2 align-items-end",
                    actionButton("export_net_png",  "Save PNG",
                                 class = "btn btn-sm btn-outline-secondary",
                                 icon  = icon("image")),
                    actionButton("export_net_jpeg", "Save JPEG",
                                 class = "btn btn-sm btn-outline-secondary",
                                 icon  = icon("image")))),
            visNetwork::visNetworkOutput("interactive_net", height = "720px"),
            div(class = "row mt-3",
                div(class = "col-md-5",
                    h5(class = "mb-2", "Top 10 Hub Nodes"),
                    DTOutput("net_hub_tbl")),
                div(class = "col-md-7",
                    uiOutput("net_edge_tbl_hdr"),
                    DTOutput("net_edge_tbl"))))
      ),
      nav_panel(
        "Network Enrichment",
        div(class = "mt-3",
            uiOutput("pathway_class_select_ui"),
            uiOutput("pathway_bubble_plot_ui"),
            div(class = "d-flex align-items-center gap-2 mt-2 mb-1",
                downloadButton("dl_net_enrich_png",  "Save plot PNG",
                               class = "btn btn-sm btn-outline-secondary", icon = icon("image")),
                downloadButton("dl_net_enrich_jpeg", "Save plot JPEG",
                               class = "btn btn-sm btn-outline-secondary", icon = icon("image")),
                downloadButton("dl_net_enrich", "Download table (.txt)",
                               class = "btn btn-sm btn-outline-primary",  icon = icon("download"))),
            div(class = "mt-2", DTOutput("tbl_pathway")))
      ),
      nav_panel(
        "Pathway Enrichment",
        div(class = "mt-3",
            uiOutput("gene_enrich_controls_ui"),
            uiOutput("gene_bubble_plot_ui"),
            div(class = "d-flex align-items-center gap-2 mt-2 mb-1",
                downloadButton("dl_pathway_enrich_png",  "Save plot PNG",
                               class = "btn btn-sm btn-outline-secondary", icon = icon("image")),
                downloadButton("dl_pathway_enrich_jpeg", "Save plot JPEG",
                               class = "btn btn-sm btn-outline-secondary", icon = icon("image")),
                downloadButton("dl_pathway_enrich", "Download table (.txt)",
                               class = "btn btn-sm btn-outline-primary",  icon = icon("download"))),
            div(class = "mt-2", DTOutput("gene_tbl_pathway")))
      ),
      nav_panel(
        "Download All",
        div(class = "mt-3",
            p("Download all result tables as a ZIP archive."),
            downloadButton("dl_all_zip", "Download results (ZIP)",
                           class = "btn btn-success"))
      )
    )
  })

  # ── CV Plot ────────────────────────────────────────────────────────────────
  output$cv_plot_or_msg <- renderUI({
    req(rv$result)   # silently return NULL before any analysis is run
    if (is.null(rv$result$CVerrors)) {
      return(div(class = "alert alert-warning",
                 "CV was not run or CVerrors not available ",
                 "(possibly because Cross-Validate was disabled)."))
    }
    tagList(
      div(class = "d-flex gap-2 mb-2",
          downloadButton("dl_cv_png",  "Save PNG",
                         class = "btn btn-sm btn-outline-secondary", icon = icon("image")),
          downloadButton("dl_cv_jpeg", "Save JPEG",
                         class = "btn btn-sm btn-outline-secondary", icon = icon("image"))),
      plotOutput("cv_plot", height = "400px", click = "cv_plot_click"),
      div(class = "d-flex align-items-end gap-3 mt-2",
          div(style = "flex: 1;",
              sliderInput("thres_slider", "Adjust threshold (λ):",
                          min = 0, max = 6, value = 2, step = 0.05,
                          width = "100%")),
          div(numericInput("thres_numeric", "λ value:",
                           value = 2, min = 0, step = 0.05, width = "100px")),
          div(class = "mb-3",
              actionButton("apply_threshold", "Update Edges & Network",
                           class = "btn btn-primary btn-sm",
                           icon  = icon("sync")))),
      uiOutput("thres_stats_ui")
    )
  })

  output$cv_plot <- renderPlot({
    req(rv$result, rv$result$CVerrors)
    cv <- rv$result$CVerrors

    # Column detection (case-insensitive)
    thres_col <- grep("threshold", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    err_col   <- grep("^CVerror$", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(err_col)) {
      err_col <- grep("cverror", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    }
    if (is.na(thres_col)) thres_col <- colnames(cv)[1L]
    if (is.na(err_col))   err_col   <- colnames(cv)[2L]

    x   <- cv[[thres_col]]
    y   <- cv[[err_col]]
    opt <- rv$result$minThres

    par(mar = c(5, 5.5, 3.5, 2), cex.axis = 1.2, cex.lab = 1.3, cex.main = 1.4)
    plot(x, y,
         type = "b", pch = 19, col = "#1a6faf", lwd = 2, cex = 1.0,
         xlab = "Threshold (λ)",
         ylab = "Mean CV Misclassification Error",
         main = "Cross-Validation: Threshold Selection",
         ylim = c(0, 1),
         bty  = "l")

    # Additional class-specific error lines
    extra_cols <- grep("^CVerror_", colnames(cv), value = TRUE)
    palette_ex <- c("#e74c3c","#2ecc71","#9b59b6","#f39c12","#1abc9c")
    for (i in seq_along(extra_cols)) {
      lines(x, cv[[extra_cols[i]]],
            col = palette_ex[(i - 1L) %% length(palette_ex) + 1L],
            lty = 2L, lwd = 1.5)
    }

    if (!is.null(opt) && is.numeric(opt) && !is.na(opt)) {
      abline(v = opt, col = "#c0392b", lty = 2L, lwd = 2L)
      legend_txt <- c("Overall CV error",
                      extra_cols,
                      paste0("Optimal λ = ", round(opt, 3L)))
      legend_col <- c("#1a6faf", palette_ex[seq_along(extra_cols)], "#c0392b")
      legend_lty <- c(1L, rep(2L, length(extra_cols)), 2L)
      sel_thres <- input$thres_slider
      if (!is.null(sel_thres) && !is.na(sel_thres)) {
        abline(v = sel_thres, col = "#27ae60", lty = 2L, lwd = 2L)
        legend_txt <- c(legend_txt, paste0("Selected λ = ", round(sel_thres, 2L)))
        legend_col <- c(legend_col, "#27ae60")
        legend_lty <- c(legend_lty, 2L)
      }
      legend("topright", legend = legend_txt, col = legend_col,
             lty = legend_lty, lwd = 2L, bty = "n", cex = 1.05)
    }
    grid(col = "grey90", lty = 1L)
  })

  # ── CV plot click: snap threshold slider to nearest x-value ────────────────
  observeEvent(input$cv_plot_click, {
    req(rv$result, rv$result$CVerrors)
    cv        <- rv$result$CVerrors
    thres_col <- grep("threshold", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(thres_col)) thres_col <- colnames(cv)[1L]
    clicked_x <- input$cv_plot_click$x
    if (is.null(clicked_x) || is.na(clicked_x)) return()
    x_vals    <- cv[[thres_col]]
    nearest   <- x_vals[which.min(abs(x_vals - clicked_x))]
    updateSliderInput(session, "thres_slider", value = nearest)
  })

  # ── Edges table — shows all edges after analysis; filtered after button click
  output$tbl_edges <- DT::renderDataTable({
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) return(DT::datatable(data.frame()))
    fe <- if (isTRUE(input$apply_threshold > 0)) filtered_edges_applied() else rv$result$EdgesSelected
    if (is.null(fe) || nrow(fe) == 0) return(NULL)
    DT::datatable(fe,
                  options   = list(pageLength = 10L, scrollX = TRUE),
                  rownames  = FALSE,
                  filter    = "top",
                  class     = "table-striped table-hover table-sm")
  })

  output$dl_edges <- downloadHandler(
    filename = function() paste0("EdgesSelected_", Sys.Date(), ".txt"),
    content  = function(file) {
      fe <- filtered_edges_applied()
      write.table(if (!is.null(fe)) fe else rv$result$EdgesSelected,
                  file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # ── CV plot download handlers ───────────────────────────────────────────────
  .render_cv_plot <- function() {
    req(rv$result, rv$result$CVerrors)
    cv        <- rv$result$CVerrors
    thres_col <- grep("threshold", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    err_col   <- grep("^CVerror$",  colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(err_col))   err_col   <- grep("cverror", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(thres_col)) thres_col <- colnames(cv)[1L]
    if (is.na(err_col))   err_col   <- colnames(cv)[2L]
    x <- cv[[thres_col]]; y <- cv[[err_col]]; opt <- rv$result$minThres
    par(mar = c(5, 5.5, 3.5, 2), cex.axis = 1.2, cex.lab = 1.3, cex.main = 1.4)
    plot(x, y, type = "b", pch = 19, col = "#1a6faf", lwd = 2, cex = 1.0,
         xlab = "Threshold (λ)", ylab = "Mean CV Misclassification Error",
         main = "Cross-Validation: Threshold Selection", ylim = c(0, 1), bty = "l")
    extra_cols <- grep("^CVerror_", colnames(cv), value = TRUE)
    palette_ex <- c("#e74c3c","#2ecc71","#9b59b6","#f39c12","#1abc9c")
    for (i in seq_along(extra_cols))
      lines(x, cv[[extra_cols[i]]], col = palette_ex[(i - 1L) %% length(palette_ex) + 1L], lty = 2L, lwd = 1.5)
    if (!is.null(opt) && is.numeric(opt) && !is.na(opt)) {
      abline(v = opt, col = "#c0392b", lty = 2L, lwd = 2L)
      sel_thres <- isolate(input$thres_slider)
      legend_txt <- c("Overall CV error", extra_cols, paste0("Optimal λ = ", round(opt, 3L)))
      legend_col <- c("#1a6faf", palette_ex[seq_along(extra_cols)], "#c0392b")
      legend_lty <- c(1L, rep(2L, length(extra_cols)), 2L)
      if (!is.null(sel_thres) && !is.na(sel_thres)) {
        abline(v = sel_thres, col = "#27ae60", lty = 2L, lwd = 2L)
        legend_txt <- c(legend_txt, paste0("Selected λ = ", round(sel_thres, 2L)))
        legend_col <- c(legend_col, "#27ae60"); legend_lty <- c(legend_lty, 2L)
      }
      legend("topright", legend = legend_txt, col = legend_col, lty = legend_lty, lwd = 2L, bty = "n", cex = 1.05)
    }
    grid(col = "grey90", lty = 1L)
  }

  output$dl_cv_png <- downloadHandler(
    filename = function() paste0("CV_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      png(file, width = 1500, height = 1000, res = 200)
      .render_cv_plot()
      dev.off()
    }
  )

  output$dl_cv_jpeg <- downloadHandler(
    filename = function() paste0("CV_plot_", Sys.Date(), ".jpg"),
    content  = function(file) {
      jpeg(file, width = 1500, height = 1000, res = 200, quality = 95)
      .render_cv_plot()
      dev.off()
    }
  )

  # ── Network PNG/JPEG export (browser-side canvas via custom JS handler) ──────
  observeEvent(input$export_net_png, {
    session$sendCustomMessage("exportVisNetwork", list(
      id         = "interactive_net",
      type       = "png",
      name       = paste0("network_", Sys.Date()),
      background = "white"
    ))
  })

  observeEvent(input$export_net_jpeg, {
    session$sendCustomMessage("exportVisNetwork", list(
      id         = "interactive_net",
      type       = "jpeg",
      name       = paste0("network_", Sys.Date()),
      background = "white"
    ))
  })

  # Store the network canvas data sent back from JS so it can be included in ZIP
  observeEvent(input$net_canvas_b64, {
    b64_str <- input$net_canvas_b64
    if (!is.character(b64_str) || !nzchar(b64_str)) return()
    b64_data <- sub("^data:[^;]+;base64,", "", b64_str)
    tryCatch({
      tmp <- file.path(tempdir(), "network_export_for_zip.png")
      writeBin(jsonlite::base64_dec(b64_data), tmp)
      rv$net_export_path <- tmp
    }, error = function(e) NULL)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # ── Sample classification download ──────────────────────────────────────────
  output$dl_sample_class <- downloadHandler(
    filename = function() paste0("SampleClassification_", Sys.Date(), ".txt"),
    content  = function(file) {
      sp <- if (!is.null(rv$result)) rv$result$SampleProbabilities else NULL
      if (is.null(sp)) sp <- data.frame()
      write.table(sp, file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # ── Network enrichment download ─────────────────────────────────────────────
  output$dl_net_enrich <- downloadHandler(
    filename = function() paste0("NetworkEnrichment_", Sys.Date(), ".txt"),
    content  = function(file) {
      df <- tryCatch(isolate(pathway_data()), error = function(e) NULL)
      if (is.null(df)) df <- data.frame()
      write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # ── Pathway enrichment download ─────────────────────────────────────────────
  output$dl_pathway_enrich <- downloadHandler(
    filename = function() paste0("PathwayEnrichment_", Sys.Date(), ".txt"),
    content  = function(file) {
      df <- tryCatch(isolate(gene_pathway_data()), error = function(e) NULL)
      if (is.null(df)) df <- data.frame()
      write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # ── Bubble plot PNG/JPEG download helpers ──────────────────────────────────
  .draw_net_enrich_plot <- function(df, pathway_class_val, pval_cutoff, left_mar = 18L) {
    if (is.null(df) || nrow(df) == 0) {
      par(mar = c(1,1,1,1)); plot.new()
      text(0.5, 0.5, "No pathways found in the enrichment results.", cex = 1.2, col = "#666")
      return(invisible())
    }
    pval_col  <- attr(df, "pval_col") %||% "HypergeoPval"
    fallback  <- isTRUE(attr(df, "fallback"))
    df <- head(df, 25L)
    df$neg_log10p <- -log10(df[[pval_col]])
    wrap_lab <- function(x) paste(strwrap(x, width = 45L), collapse = "\n")
    df$label  <- sapply(as.character(df$PathAnnotation), wrap_lab, USE.NAMES = FALSE)
    df$label  <- factor(df$label, levels = rev(unique(df$label)))
    df$col    <- ifelse(df$Direction == "up", "#e74c3c", "#3498db")
    df$size_val <- pmax(1, df$Enriched_Genesize %||% df$Enriched_Edgesize %||% 1)
    xlab_expr  <- if (pval_col == "AdjPval") expression(-log[10](adj~p)) else expression(-log[10](p))
    main_title <- paste0(pathway_class_val, " — Network Enrichment")
    if (fallback) main_title <- paste0(main_title, "\n(no significant pathways — top 15 shown)")
    all_lines   <- unlist(strsplit(as.character(levels(df$label)), "\n"))
    max_lc      <- max(nchar(all_lines), na.rm = TRUE)
    line_counts <- sapply(as.character(levels(df$label)), function(x) length(strsplit(x, "\n")[[1L]]))
    y_step      <- max(line_counts)
    y_positions <- cumsum(c(1, head(pmax(line_counts, 1.5), -1)))
    df$y_pos    <- y_positions[as.integer(df$label)]
    par(mar = c(5, left_mar, 3, 2),
        cex.axis = 1.1, cex.lab = 1.1, cex.main = 1.4, lheight = 0.8)
    plot(df$neg_log10p, df$y_pos, pch = 21, bg = df$col, col = "#333333",
         cex = 1.2 + sqrt(df$size_val) * 0.55,
         xlab = xlab_expr, ylab = "", yaxt = "n", bty = "l",
         xlim = c(0, max(df$neg_log10p) * 1.1),
         ylim = c(0.5, max(y_positions) + max(y_step, 1.5) * 0.5),
         main = main_title)
    axis(2, at = y_positions, labels = levels(df$label), las = 1, cex.axis = 1.1, gap.axis = 0)
    abline(v = -log10(pval_cutoff), col = "#888", lty = 2)
    legend("bottomright", legend = c("Up-regulated", "Down-regulated"),
           pt.bg = c("#e74c3c", "#3498db"), pch = 21, pt.cex = 1.6,
           col = "#333333", bty = "n", cex = 1.1)
  }

  .draw_gene_enrich_plot <- function(df, group_val, pval_cutoff, left_mar = 18L) {
    if (is.null(df) || nrow(df) == 0) {
      par(mar = c(1,1,1,1)); plot.new()
      text(0.5, 0.5,
           "No pathway enrichment results.\nEnsure analysis is run with Data X genes and a pathway file is loaded.",
           cex = 1.1, col = "#666", adj = 0.5)
      return(invisible())
    }
    pval_col  <- attr(df, "pval_col") %||% "HypergeoPval"
    fallback  <- isTRUE(attr(df, "fallback"))
    df <- head(df, 25L)
    df$neg_log10p <- -log10(df[[pval_col]])
    wrap_lab <- function(x) paste(strwrap(x, width = 45L), collapse = "\n")
    df$label  <- sapply(as.character(df$PathAnnotation), wrap_lab, USE.NAMES = FALSE)
    df$label  <- factor(df$label, levels = rev(unique(df$label)))
    df$col    <- ifelse(df$Direction == "up", "#e74c3c", "#3498db")
    df$size_val <- pmax(1, df$Enriched_Genesize)
    xlab_expr  <- if (pval_col == "AdjPval") expression(-log[10](adj~p)) else expression(-log[10](p))
    main_title <- paste0(group_val, " — Gene Pathway Enrichment")
    if (fallback) main_title <- paste0(main_title, "\n(no significant pathways — top 15 shown)")
    all_lines   <- unlist(strsplit(as.character(levels(df$label)), "\n"))
    max_lc      <- max(nchar(all_lines), na.rm = TRUE)
    line_counts <- sapply(as.character(levels(df$label)), function(x) length(strsplit(x, "\n")[[1L]]))
    y_step      <- max(line_counts)
    y_positions <- cumsum(c(1, head(pmax(line_counts, 1.5), -1)))
    df$y_pos    <- y_positions[as.integer(df$label)]
    par(mar = c(5, left_mar, 3, 2),
        cex.axis = 1.1, cex.lab = 1.1, cex.main = 1.4, lheight = 0.8)
    plot(df$neg_log10p, df$y_pos, pch = 21, bg = df$col, col = "#333333",
         cex = 1.2 + sqrt(df$size_val) * 0.55,
         xlab = xlab_expr, ylab = "", yaxt = "n", bty = "l",
         xlim = c(0, max(df$neg_log10p) * 1.1),
         ylim = c(0.5, max(y_positions) + max(y_step, 1.5) * 0.5),
         main = main_title)
    axis(2, at = y_positions, labels = levels(df$label), las = 1, cex.axis = 1.1, gap.axis = 0)
    abline(v = -log10(pval_cutoff), col = "#888", lty = 2)
    legend("bottomright", legend = c("Up-regulated", "Down-regulated"),
           pt.bg = c("#e74c3c", "#3498db"), pch = 21, pt.cex = 1.6,
           col = "#333333", bty = "n", cex = 1.1)
  }

  output$dl_net_enrich_png <- downloadHandler(
    filename = function() paste0("NetworkEnrichment_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      df   <- tryCatch(isolate(pathway_data()),   error = function(e) NULL)
      cls  <- isolate(input$pathway_class) %||% ""
      pcut <- isolate(input$pathway_pval)  %||% 0.05
      dims <- .bubble_dims(df)
      png(file, width = dims$width, height = dims$height, res = 200)
      .draw_net_enrich_plot(df, cls, pcut, left_mar = dims$left_mar)
      dev.off()
    }
  )
  output$dl_net_enrich_jpeg <- downloadHandler(
    filename = function() paste0("NetworkEnrichment_plot_", Sys.Date(), ".jpg"),
    content  = function(file) {
      df   <- tryCatch(isolate(pathway_data()),   error = function(e) NULL)
      cls  <- isolate(input$pathway_class) %||% ""
      pcut <- isolate(input$pathway_pval)  %||% 0.05
      dims <- .bubble_dims(df)
      jpeg(file, width = dims$width, height = dims$height, res = 200, quality = 95)
      .draw_net_enrich_plot(df, cls, pcut, left_mar = dims$left_mar)
      dev.off()
    }
  )
  output$dl_pathway_enrich_png <- downloadHandler(
    filename = function() paste0("PathwayEnrichment_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      df   <- tryCatch(isolate(gene_pathway_data()), error = function(e) NULL)
      grp  <- isolate(input$gene_enrich_group) %||% ""
      pcut <- isolate(input$gene_enrich_pval)  %||% 0.05
      dims <- .bubble_dims(df)
      png(file, width = dims$width, height = dims$height, res = 200)
      .draw_gene_enrich_plot(df, grp, pcut, left_mar = dims$left_mar)
      dev.off()
    }
  )
  output$dl_pathway_enrich_jpeg <- downloadHandler(
    filename = function() paste0("PathwayEnrichment_plot_", Sys.Date(), ".jpg"),
    content  = function(file) {
      df   <- tryCatch(isolate(gene_pathway_data()), error = function(e) NULL)
      grp  <- isolate(input$gene_enrich_group) %||% ""
      pcut <- isolate(input$gene_enrich_pval)  %||% 0.05
      dims <- .bubble_dims(df)
      jpeg(file, width = dims$width, height = dims$height, res = 200, quality = 95)
      .draw_gene_enrich_plot(df, grp, pcut, left_mar = dims$left_mar)
      dev.off()
    }
  )

  # ── Sample probabilities + confusion matrix ────────────────────────────────
  output$tbl_probs <- DT::renderDataTable({
    if (is.null(rv$result) || is.null(rv$result$SampleProbabilities)) return(DT::datatable(data.frame()))
    DT::datatable(rv$result$SampleProbabilities,
                  options   = list(pageLength = 10L, scrollX = TRUE),
                  rownames  = FALSE,
                  filter    = "top",
                  class     = "table-striped table-hover table-sm")
  })

  output$confusion_matrix_ui <- renderUI({
    if (is.null(rv$result)) return(NULL)
    sp <- rv$result$SampleProbabilities
    if (is.null(sp)) return(NULL)

    tc <- grep("TrueClass|True_Class",      colnames(sp), ignore.case = TRUE, value = TRUE)[1L]
    pc <- grep("PredictedClass|Pred_Class", colnames(sp), ignore.case = TRUE, value = TRUE)[1L]

    if (is.na(tc) || is.na(pc)) {
      return(div(class = "alert alert-warning",
                 "Could not find TrueClass / PredictedClass columns for confusion matrix."))
    }

    conf <- table(True = sp[[tc]], Predicted = sp[[pc]])
    acc  <- round(100 * sum(diag(conf)) / sum(conf), 1L)

    tagList(
      div(class = "alert alert-info d-inline-block px-3 py-2",
          tags$b(paste0("Overall Accuracy: ", acc, "%"))),
      div(
        class = "table-responsive mt-2",
        tags$table(
          class = "table table-bordered table-sm w-auto",
          tags$thead(
            tags$tr(
              tags$th("True ↓ Predicted →"),
              lapply(colnames(conf), function(cn) tags$th(cn, class = "text-center"))
            )
          ),
          tags$tbody(
            lapply(rownames(conf), function(rn) {
              tags$tr(
                tags$th(rn, scope = "row"),
                lapply(colnames(conf), function(cn) {
                  val <- conf[rn, cn]
                  cls <- if (rn == cn && val > 0) "table-success"
                  else if (val > 0) "table-danger"
                  else ""
                  tags$td(val, class = cls, style = "text-align:center;")
                })
              )
            })
          )
        )
      )
    )
  })

  # ── filtered_edges: live reactive for stats card / CV plot ─────────────────
  filtered_edges <- reactive({
    req(rv$result, rv$result$EdgesSelected)
    es        <- rv$result$EdgesSelected
    es_cn     <- names(es); if (is.null(es_cn)) es_cn <- colnames(es)
    abs_cols  <- grep("_absdiknew$", es_cn, value = TRUE, ignore.case = TRUE)
    if (length(abs_cols) == 0) return(es)
    opt       <- if (!is.null(rv$result$minThres)) rv$result$minThres else 0
    max_orig  <- apply(es[, abs_cols, drop = FALSE], 1, max, na.rm = TRUE) + opt
    thres_val <- if (!is.null(input$thres_slider)) input$thres_slider else opt
    es[max_orig > thres_val, , drop = FALSE]
  })

  # ── filtered_edges_applied: only updates on button click ───────────────────
  # Delegates to filtered_edges() to avoid re-running the apply() computation.
  filtered_edges_applied <- eventReactive(
    input$apply_threshold,
    { req(filtered_edges()); filtered_edges() },
    ignoreNULL = TRUE
  )

  # ── Update slider bounds/value when new result arrives ─────────────────────
  observeEvent(rv$result, {
    req(rv$result, rv$result$EdgesSelected)
    es       <- rv$result$EdgesSelected
    cn       <- names(es)
    if (is.null(cn)) cn <- colnames(es)
    abs_cols <- grep("_absdiknew$", cn, value = TRUE, ignore.case = TRUE)
    # Store class names for the edge-colour dropdown
    rv$edge_classes <- sub("(?i)_absdiknew$", "", abs_cols, perl = TRUE)
    opt      <- if (!is.null(rv$result$minThres)) rv$result$minThres else 2
    if (length(abs_cols) > 0) {
      max_orig <- apply(es[, abs_cols, drop = FALSE], 1, max, na.rm = TRUE) + opt
      new_max  <- ceiling(max(max_orig, na.rm = TRUE) * 10) / 10
    } else {
      new_max <- 6
    }
    updateSliderInput(session, "thres_slider",
                      min   = 0,
                      max   = new_max,
                      value = opt)
  })

  # ── Threshold stats info card ───────────────────────────────────────────────
  output$thres_stats_ui <- renderUI({
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) return(NULL)
    req(input$thres_slider)
    n_edges <- nrow(filtered_edges())
    cv_err_pct <- NA_real_
    if (!is.null(rv$result$CVerrors)) {
      cv <- rv$result$CVerrors
      tc <- grep("threshold", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
      ec <- grep("^CVerror$",  colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
      if (is.na(ec)) ec <- grep("cverror", colnames(cv), ignore.case = TRUE, value = TRUE)[1L]
      if (!is.na(tc) && !is.na(ec)) {
        closest    <- which.min(abs(cv[[tc]] - input$thres_slider))
        cv_err_pct <- round(cv[[ec]][closest] * 100, 1)
      }
    }
    div(
      class = "alert alert-info mt-2",
      style = "display:inline-block; padding:6px 12px;",
      HTML(sprintf(
        "At λ = %.2f: <strong>%d edges</strong> selected &middot; CV error = %s",
        input$thres_slider,
        n_edges,
        if (is.na(cv_err_pct)) "N/A" else paste0(cv_err_pct, "%")
      ))
    )
  })

  # ── Dynamic badge for Selected Edges tab ──────────────────────────────────
  output$edges_count_badge <- renderText({
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) return(NULL)
    n <- if (isTRUE(input$apply_threshold > 0) && !is.null(filtered_edges_applied())) {
      nrow(filtered_edges_applied())
    } else {
      nrow(rv$result$EdgesSelected)
    }
    paste0("(", n, " edges)")
  })

  # ── Bidirectional sync: slider ↔ numeric input ─────────────────────────────
  observeEvent(input$thres_slider, {
    req(input$thres_slider)
    if (!isTRUE(all.equal(as.numeric(input$thres_numeric %||% -1),
                          as.numeric(input$thres_slider)))) {
      updateNumericInput(session, "thres_numeric", value = input$thres_slider)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$thres_numeric, {
    req(input$thres_numeric, is.numeric(input$thres_numeric), !is.na(input$thres_numeric))
    if (!isTRUE(all.equal(as.numeric(input$thres_slider %||% -1),
                          as.numeric(input$thres_numeric)))) {
      updateSliderInput(session, "thres_slider", value = input$thres_numeric)
    }
  }, ignoreInit = TRUE)

  # ── Network group coloring selector ────────────────────────────────────────
  output$net_group_select_ui <- renderUI({
    if (is.null(rv$result)) return(NULL)
    classes <- rv$edge_classes   # populated in observeEvent(rv$result)
    if (length(classes) == 0) return(NULL)
    choices <- c("Edge magnitude" = "magnitude",
                 setNames(classes, paste("Direction:", classes)))
    div(class = "mb-2",
        selectInput("net_group_color", "Color edges by:",
                    choices = choices, width = "220px"))
  })

  # ── Network enrichment controls ────────────────────────────────────────────
  output$pathway_class_select_ui <- renderUI({
    # Use R-computed results if available, else fall back to C++ results
    enrich_src <- if (!is.null(rv$net_enrich_result)) rv$net_enrich_result
                  else if (!is.null(rv$result))        rv$result
                  else NULL
    req(!is.null(enrich_src))
    eu      <- enrich_src$Enrichment_up
    groups  <- if (!is.null(eu) && length(eu) > 0) names(eu) else character(0)
    if (length(groups) == 0)
      return(div(class = "alert alert-warning",
                 "No enrichment results found. Use 'Run Network Enrichment' in the left panel after loading a pathway file."))
    div(class = "d-flex flex-wrap gap-3 align-items-center mb-2",
        selectInput("pathway_class", "Group:", choices = groups, width = "150px"),
        radioButtons("pathway_dir", "Direction:",
                     choices = c("Up-regulated" = "up", "Down-regulated" = "down", "Both" = "both"),
                     inline = TRUE),
        numericInput("pathway_pval", "p-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01, width = "120px"),
        numericInput("pathway_fdr",  "FDR cutoff:",     value = 0.25, min = 0, max = 1, step = 0.01, width = "120px"),
        selectInput("net_enrich_correct", "Correction:",
                    choices = c("None" = "none", "BH (FDR)" = "BH",
                                "Bonferroni" = "bonferroni", "q-value" = "qvalue"),
                    selected = "BH", width = "140px"))
  })

  # ── Network enrichment bubble plot data ────────────────────────────────────
  pathway_data <- reactive({
    enrich_src <- if (!is.null(rv$net_enrich_result)) rv$net_enrich_result
                  else if (!is.null(rv$result))        rv$result
                  else NULL
    req(!is.null(enrich_src), input$pathway_class)
    cls     <- input$pathway_class
    dir     <- if (!is.null(input$pathway_dir)) input$pathway_dir else "both"
    pcut    <- if (!is.null(input$pathway_pval) && !is.na(input$pathway_pval)) input$pathway_pval else 0.05
    fdr_cut <- if (!is.null(input$pathway_fdr)  && !is.na(input$pathway_fdr))  input$pathway_fdr  else 0.25
    correct <- input$net_enrich_correct %||% "BH"

    get_enrich <- function(direction) {
      lst <- if (direction == "up") enrich_src$Enrichment_up else enrich_src$Enrichment_down
      df  <- tryCatch(lst[[cls]], error = function(e) NULL)
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
      df$Direction <- direction
      # Apply correction if not already present (C++ results lack AdjPval)
      if (!"AdjPval" %in% names(df)) {
        pv <- df$HypergeoPval
        valid <- !is.na(pv)
        p_adj <- pv
        if (correct != "none" && sum(valid) > 1L) {
          p_adj[valid] <- if (correct == "qvalue") {
            tryCatch(qvalue::qvalue(pv[valid], pi0.method = "smoother")$qvalues,
                     error = function(e) p.adjust(pv[valid], method = "BH"))
          } else {
            p.adjust(pv[valid], method = correct)
          }
        }
        df$AdjPval <- p_adj
      }
      df
    }

    dfs <- if (dir == "both") {
      Filter(Negate(is.null), list(get_enrich("up"), get_enrich("down")))
    } else {
      Filter(Negate(is.null), list(get_enrich(dir)))
    }
    if (length(dfs) == 0) return(NULL)
    df <- do.call(rbind, dfs)

    pval_col <- if (correct != "none" && "AdjPval" %in% names(df)) "AdjPval" else "HypergeoPval"
    df <- df[!is.na(df$HypergeoPval), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[order(df[[pval_col]]), , drop = FALSE]
    attr(df, "pval_col") <- pval_col
    keep <- df$HypergeoPval <= pcut
    if (correct != "none" && "AdjPval" %in% names(df)) keep <- keep & (df$AdjPval <= fdr_cut)
    sig <- df[keep, , drop = FALSE]
    if (nrow(sig) > 0) {
      attr(sig, "pval_col") <- pval_col
      return(sig)
    }
    out <- head(df, 15L)
    attr(out, "fallback")  <- TRUE
    attr(out, "pval_col")  <- pval_col
    out
  })

  # Helper: compute bubble plot height for interactive display
  .bubble_height <- function(df, wrap_width = 45L) {
    if (is.null(df) || nrow(df) == 0L) return(480L)
    wrap_lab <- function(x) paste(strwrap(x, width = wrap_width), collapse = "\n")
    labs <- sapply(as.character(head(df$PathAnnotation, 25L)), wrap_lab, USE.NAMES = FALSE)
    y_step  <- max(sapply(labs, function(x) length(strsplit(x, "\n")[[1L]])), 1L)
    n_items <- length(unique(labs))
    px_per_item <- 10L * y_step + 5L
    max(480L, as.integer((n_items + 1L) * px_per_item + 120L))
  }

  # Helper: compute PNG width/height/left-margin for saved bubble plots (res=200)
  # At res=200, cex.axis=1.1, lheight=0.8: 1 margin line ≈ 27px.
  # We budget 1500px for the plot body, and set width = left_px + 1500.
  .bubble_dims <- function(df, wrap_width = 45L) {
    if (is.null(df) || nrow(df) == 0L)
      return(list(width = 2000L, height = 1000L, left_mar = 18L))
    wrap_lab <- function(x) paste(strwrap(x, width = wrap_width), collapse = "\n")
    labs <- sapply(as.character(head(df$PathAnnotation, 25L)), wrap_lab, USE.NAMES = FALSE)
    all_lines   <- unlist(strsplit(labs, "\n"))
    max_lc      <- max(nchar(all_lines), na.rm = TRUE)
    line_counts <- sapply(labs, function(x) length(strsplit(x, "\n")[[1L]]))
    y_step      <- max(line_counts, 1L)
    n_items     <- length(unique(labs))
    left_mar    <- max(14L, as.integer(max_lc * 0.67 + 2L))
    left_px     <- left_mar * 27L
    width_px    <- max(1700L, left_px + 1500L)
    px_per_item <- 23L * y_step + 10L
    height_px   <- max(1000L, as.integer((n_items + 1L) * px_per_item + 250L))
    list(width = width_px, height = height_px, left_mar = left_mar)
  }

  output$pathway_bubble_plot_ui <- renderUI({
    h <- tryCatch(.bubble_height(pathway_data()), error = function(e) 480L)
    div(style = "overflow-x: auto;",
        plotOutput("pathway_bubble_plot", height = paste0(h, "px"), width = "1000px"))
  })

  output$gene_bubble_plot_ui <- renderUI({
    h <- tryCatch(.bubble_height(gene_pathway_data()), error = function(e) 480L)
    div(style = "overflow-x: auto;",
        plotOutput("gene_bubble_plot", height = paste0(h, "px"), width = "1000px"))
  })

  output$pathway_bubble_plot <- renderPlot({
    df <- pathway_data()
    if (is.null(df) || nrow(df) == 0) {
      par(mar = c(1, 1, 1, 1))
      plot.new()
      text(0.5, 0.5, "No pathways found in the enrichment results.",
           cex = 1.2, col = "#666")
      return()
    }

    pval_col <- attr(df, "pval_col") %||% "HypergeoPval"
    fallback <- isTRUE(attr(df, "fallback"))

    df <- head(df, 25L)
    df$neg_log10p <- -log10(df[[pval_col]])
    wrap_lab <- function(x) paste(strwrap(x, width = 45L), collapse = "\n")
    df$label <- sapply(as.character(df$PathAnnotation), wrap_lab, USE.NAMES = FALSE)
    df$label <- factor(df$label, levels = rev(unique(df$label)))
    df$col   <- ifelse(df$Direction == "up", "#e74c3c", "#3498db")
    df$size_val <- pmax(1, df$Enriched_Genesize %||% df$Enriched_Edgesize %||% 1)

    n_paths    <- nlevels(df$label)
    xlab_expr  <- if (pval_col == "AdjPval") expression(-log[10](adj~p)) else expression(-log[10](p))
    main_title <- paste0(input$pathway_class, " — Network Enrichment")
    if (fallback) main_title <- paste0(main_title, "\n(no significant pathways — top 15 shown)")

    all_lines   <- unlist(strsplit(as.character(levels(df$label)), "\n"))
    max_lc      <- max(nchar(all_lines), na.rm = TRUE)
    line_counts <- sapply(as.character(levels(df$label)), function(x) length(strsplit(x, "\n")[[1L]]))
    y_step      <- max(line_counts)
    y_positions <- cumsum(c(1, head(pmax(line_counts, 1.5), -1)))
    df$y_pos    <- y_positions[as.integer(df$label)]

    par(mar = c(5, max(18, max_lc * 0.61 + 2), 3, 2),
        cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.56, lheight = 0.8)
    plot(df$neg_log10p, df$y_pos,
         pch = 21, bg = df$col, col = "#333333",
         cex = 1.2 + sqrt(df$size_val) * 0.55,
         xlab = xlab_expr, ylab = "",
         yaxt = "n", bty = "l",
         xlim = c(0, max(df$neg_log10p) * 1.1),
         ylim = c(0.5, max(y_positions) + max(y_step, 1.5) * 0.5),
         main = main_title)
    axis(2, at = y_positions, labels = levels(df$label), las = 1, cex.axis = 1.3, gap.axis = 0)
    abline(v = -log10(input$pathway_pval %||% 0.05), col = "#888", lty = 2)
    legend("bottomright", legend = c("Up-regulated", "Down-regulated"),
           pt.bg = c("#e74c3c", "#3498db"), pch = 21, pt.cex = 1.8,
           col = "#333333", bty = "n", cex = 1.3)
  })

  output$tbl_pathway <- DT::renderDataTable({
    df <- pathway_data()
    if (is.null(df)) return(DT::datatable(data.frame()))
    fallback <- isTRUE(attr(df, "fallback"))
    show_cols <- intersect(c("Pathway","PathAnnotation","Direction","HypergeoPval","AdjPval",
                              "Enriched_Edgesize","NumGene_inPathway","Geneset_Enriched"),
                           colnames(df))
    caption <- if (fallback)
      "No significant pathways at the selected cutoff — showing top 15 by p-value."
    else NULL
    DT::datatable(df[, show_cols, drop = FALSE],
                  caption  = caption,
                  options  = list(pageLength = 15L, scrollX = TRUE),
                  rownames = FALSE, filter = "top",
                  class    = "table-striped table-hover table-sm")
  })

  # ── Gene-level pathway enrichment controls ────────────────────────────────
  output$gene_enrich_controls_ui <- renderUI({
    req(rv$result, rv$result$EdgesSelected)   # silently return NULL before analysis
    if (is.null(rv$pathway))
      return(div(class = "alert alert-warning",
                 "Load a pathway annotation file to use gene-level pathway enrichment."))
    es     <- rv$result$EdgesSelected
    groups <- sub("(?i)_absdiknew$", "",
                  grep("_absdiknew$", names(es), value = TRUE, ignore.case = TRUE),
                  perl = TRUE)
    if (length(groups) == 0)
      return(div(class = "alert alert-warning", "No group columns found in selected edges."))
    pw_classes <- sort(unique(sub(":.*", "", as.character(rv$pathway[[2L]]))))
    pw_choices  <- c("All classes" = "ALL", setNames(pw_classes, pw_classes))
    has_Y       <- !is.null(rv$data_Y)
    data_choices <- if (has_Y)
      c("Data X" = "X", "Data Y" = "Y", "Both" = "both")
    else
      c("Data X" = "X")
    div(class = "d-flex flex-wrap gap-3 align-items-center mb-2",
        selectInput("gene_enrich_group",   "Group:",    choices = groups,   width = "150px"),
        radioButtons("gene_enrich_data",   "Include data:",
                     choices = data_choices, inline = TRUE),
        radioButtons("gene_enrich_dir",    "Direction:",
                     choices = c("Up-regulated" = "up", "Down-regulated" = "down", "Both" = "both"),
                     inline = TRUE),
        selectInput("gene_enrich_class",   "Pathway class:", choices = pw_choices, width = "200px"),
        numericInput("gene_enrich_pval",   "p-value cutoff:",
                     value = 0.05, min = 0, max = 1, step = 0.01, width = "130px"),
        numericInput("gene_enrich_fdr",    "FDR cutoff:",
                     value = 0.25, min = 0, max = 1, step = 0.01, width = "130px"),
        numericInput("gene_enrich_minbg",  "Min genes in pathway:",
                     value = 5L,   min = 1L, max = 200L, step = 1L, width = "160px"),
        selectInput("gene_enrich_correct", "Correction:",
                    choices = c("None" = "none", "BH (FDR)" = "BH",
                                "Bonferroni" = "bonferroni", "q-value" = "qvalue"),
                    selected = "BH", width = "140px"))
  })

  # ── Gene-level hypergeometric enrichment ──────────────────────────────────
  gene_pathway_data <- reactive({
    rv$gene_enrich_counter  # invalidate when "Run Pathway Enrichment" button is clicked
    req(rv$result, rv$pathway, input$gene_enrich_group)
    grp     <- input$gene_enrich_group
    dir     <- input$gene_enrich_dir   %||% "both"
    pcut    <- if (!is.null(input$gene_enrich_pval)  && !is.na(input$gene_enrich_pval))  input$gene_enrich_pval  else 0.05
    fdr_cut <- if (!is.null(input$gene_enrich_fdr)   && !is.na(input$gene_enrich_fdr))   input$gene_enrich_fdr   else 0.25
    minbg   <- if (!is.null(input$gene_enrich_minbg) && !is.na(input$gene_enrich_minbg)) input$gene_enrich_minbg else 5L
    cls     <- input$gene_enrich_class   %||% "ALL"
    correct <- input$gene_enrich_correct %||% "BH"

    data_sel <- input$gene_enrich_data %||% "X"
    bg_src   <- input$gene_enrich_bg   %||% "data"

    # Background: data features / all pathway genes / HGNC / custom upload
    bg_genes <- {
      if (bg_src == "custom" && (is.null(rv$bg_custom) || length(rv$bg_custom) == 0L)) {
        showNotification("No custom background loaded — using data features instead.",
                         type = "warning", duration = 5L)
        bg_src <- "data"
      }
      if (bg_src == "hgnc_all") {
        if (length(.hgnc_all) == 0L) {
          showNotification("HGNC gene list not found in package — using data features.",
                           type = "warning", duration = 5L)
          bg_X <- if (!is.null(rv$data_X)) unique(as.character(rv$data_X[[1L]])) else character(0)
          bg_X
        } else .hgnc_all
      } else if (bg_src == "hgnc_protein_product") {
        if (length(.hgnc_protein_product) == 0L) {
          showNotification("HGNC protein-product list not found — using data features.",
                           type = "warning", duration = 5L)
          bg_X <- if (!is.null(rv$data_X)) unique(as.character(rv$data_X[[1L]])) else character(0)
          bg_X
        } else .hgnc_protein_product
      } else if (bg_src == "custom") {
        rv$bg_custom
      } else if (bg_src == "pathway") {
        if (!is.null(rv$pathway)) unique(as.character(rv$pathway[[1L]])) else character(0)
      } else {
        # "data" — features from the selected data layer(s)
        bg_X <- if (!is.null(rv$data_X)) unique(as.character(rv$data_X[[1L]])) else character(0)
        bg_Y <- if (!is.null(rv$data_Y)) unique(as.character(rv$data_Y[[1L]])) else character(0)
        switch(data_sel, "X" = bg_X, "Y" = bg_Y, "both" = unique(c(bg_X, bg_Y)), bg_X)
      }
    }
    if (length(bg_genes) == 0) return(NULL)

    es <- rv$result$EdgesSelected
    if (is.null(es)) return(NULL)

    # Extract gene names from the selected data layer(s) for a given direction
    get_signal <- function(direction) {
      dir_col <- paste0(grp, "_dir")
      sub_es  <- if (direction == "both" || !dir_col %in% names(es)) es
                 else es[!is.na(es[[dir_col]]) & es[[dir_col]] == direction, , drop = FALSE]
      pat <- switch(data_sel,
        "X"    = "\\.\\.\\.datX$",
        "Y"    = "\\.\\.\\.datY$",
        "both" = "\\.\\.\\.(datX|datY)$",
        "\\.\\.\\.datX$"
      )
      nodeA <- sub("\\.\\.\\.(datX|datY)$", "", sub_es$NodeA[grepl(pat, sub_es$NodeA)])
      nodeB <- sub("\\.\\.\\.(datX|datY)$", "", sub_es$NodeB[grepl(pat, sub_es$NodeB)])
      unique(c(nodeA, nodeB))
    }

    run_hypergeo <- function(direction) {
      signal <- get_signal(direction)
      if (length(signal) == 0) return(NULL)

      pw <- rv$pathway
      if (cls != "ALL") pw <- pw[sub(":.*", "", as.character(pw[[2L]])) == cls, , drop = FALSE]
      if (nrow(pw) == 0) return(NULL)

      pw_genes_col <- as.character(pw[[1L]])
      pw_ids_col   <- as.character(pw[[2L]])
      pw_name_col  <- as.character(pw[[3L]])

      # Per-pathway gene sets (deduplicated)
      pw_list <- split(pw_genes_col, pw_ids_col)
      pw_name <- tapply(pw_name_col, pw_ids_col, function(x) x[[1L]])

      bg_set  <- bg_genes
      sig_set <- signal
      N <- length(bg_set)
      n <- length(sig_set)

      # Background pathway sizes; drop pathways with too few bg genes
      K_v <- vapply(pw_list, function(g) length(intersect(unique(g), bg_set)), integer(1L))
      keep    <- K_v >= minbg
      pw_list <- pw_list[keep]; pw_name <- pw_name[keep]; K_v <- K_v[keep]
      if (length(pw_list) == 0) return(NULL)

      k_v <- vapply(pw_list, function(g) length(intersect(unique(g), sig_set)), integer(1L))
      p_v <- phyper(k_v - 1L, K_v, N - K_v, n, lower.tail = FALSE)

      # Multiple-testing correction (per-direction, before combining)
      p_adj <- if (correct == "none" || length(p_v) <= 1L) {
        p_v
      } else if (correct == "qvalue") {
        tryCatch(
          qvalue::qvalue(p_v, pi0.method = "smoother")$qvalues,
          error = function(e) p.adjust(p_v, method = "BH")
        )
      } else {
        p.adjust(p_v, method = correct)
      }

      # Gene lists for hits
      gene_str <- vapply(pw_list, function(g) {
        hits <- intersect(unique(g), sig_set)
        if (length(hits) == 0) "" else paste(sort(hits), collapse = ";")
      }, character(1L))

      data.frame(
        Pathway           = names(pw_list),
        PathAnnotation    = as.character(pw_name),
        Direction         = direction,
        HypergeoPval      = p_v,
        AdjPval           = p_adj,
        Enriched_Genesize = k_v,
        NumGene_inPathway = K_v,
        Geneset_Enriched  = gene_str,
        stringsAsFactors  = FALSE
      )
    }

    dfs <- if (dir == "both") {
      Filter(Negate(is.null), list(run_hypergeo("up"), run_hypergeo("down")))
    } else {
      Filter(Negate(is.null), list(run_hypergeo(dir)))
    }
    if (length(dfs) == 0) return(NULL)
    df <- do.call(rbind, dfs)

    pval_col <- if (correct != "none") "AdjPval" else "HypergeoPval"
    df <- df[!is.na(df$HypergeoPval), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[order(df[[pval_col]]), , drop = FALSE]
    attr(df, "pval_col") <- pval_col
    keep <- df$HypergeoPval <= pcut
    if (correct != "none" && "AdjPval" %in% names(df)) keep <- keep & (df$AdjPval <= fdr_cut)
    sig <- df[keep, , drop = FALSE]
    if (nrow(sig) > 0) {
      attr(sig, "pval_col") <- pval_col
      return(sig)
    }
    out <- head(df, 15L)
    attr(out, "fallback")  <- TRUE
    attr(out, "pval_col")  <- pval_col
    out
  })

  output$gene_bubble_plot <- renderPlot({
    df <- gene_pathway_data()
    if (is.null(df) || nrow(df) == 0) {
      par(mar = c(1, 1, 1, 1)); plot.new()
      text(0.5, 0.5,
           "No pathway enrichment results.\nEnsure analysis is run with Data X genes and a pathway file is loaded.",
           cex = 1.1, col = "#666", adj = 0.5)
      return()
    }
    pval_col  <- attr(df, "pval_col") %||% "HypergeoPval"
    fallback  <- isTRUE(attr(df, "fallback"))
    df <- head(df, 25L)
    df$neg_log10p <- -log10(df[[pval_col]])
    wrap_lab <- function(x) paste(strwrap(x, width = 45L), collapse = "\n")
    df$label <- sapply(as.character(df$PathAnnotation), wrap_lab, USE.NAMES = FALSE)
    df$label <- factor(df$label, levels = rev(unique(df$label)))
    df$col   <- ifelse(df$Direction == "up", "#e74c3c", "#3498db")
    df$size_val <- pmax(1, df$Enriched_Genesize)
    n_paths    <- nlevels(df$label)
    xlab_expr  <- if (pval_col == "AdjPval") expression(-log[10](adj~p)) else expression(-log[10](p))
    grp_label  <- if (!is.null(input$gene_enrich_group)) input$gene_enrich_group else ""
    main_title <- paste0(grp_label, " — Gene Pathway Enrichment")
    if (fallback) main_title <- paste0(main_title, "\n(no significant pathways — top 15 shown)")

    all_lines   <- unlist(strsplit(as.character(levels(df$label)), "\n"))
    max_lc      <- max(nchar(all_lines), na.rm = TRUE)
    line_counts <- sapply(as.character(levels(df$label)), function(x) length(strsplit(x, "\n")[[1L]]))
    y_step      <- max(line_counts)
    y_positions <- cumsum(c(1, head(pmax(line_counts, 1.5), -1)))
    df$y_pos    <- y_positions[as.integer(df$label)]

    par(mar = c(5, max(18, max_lc * 0.61 + 2), 3, 2),
        cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.56, lheight = 0.8)
    plot(df$neg_log10p, df$y_pos,
         pch = 21, bg = df$col, col = "#333333",
         cex = 1.2 + sqrt(df$size_val) * 0.55,
         xlab = xlab_expr, ylab = "",
         yaxt = "n", bty = "l",
         xlim = c(0, max(df$neg_log10p) * 1.1),
         ylim = c(0.5, max(y_positions) + max(y_step, 1.5) * 0.5),
         main = main_title)
    axis(2, at = y_positions, labels = levels(df$label), las = 1, cex.axis = 1.3, gap.axis = 0)
    abline(v = -log10(input$gene_enrich_pval %||% 0.05), col = "#888", lty = 2)
    legend("bottomright", legend = c("Up-regulated", "Down-regulated"),
           pt.bg = c("#e74c3c", "#3498db"), pch = 21, pt.cex = 1.8,
           col = "#333333", bty = "n", cex = 1.3)
  })

  output$gene_tbl_pathway <- DT::renderDataTable({
    df <- gene_pathway_data()
    if (is.null(df)) return(DT::datatable(data.frame()))
    fallback  <- isTRUE(attr(df, "fallback"))
    show_cols <- intersect(c("Pathway","PathAnnotation","Direction","HypergeoPval","AdjPval",
                              "Enriched_Genesize","NumGene_inPathway","Geneset_Enriched"),
                           colnames(df))
    caption <- if (fallback) "No significant pathways at the selected cutoff — showing top 15 by p-value." else NULL
    DT::datatable(df[, show_cols, drop = FALSE],
                  caption  = caption,
                  options  = list(pageLength = 15L, scrollX = TRUE),
                  rownames = FALSE, filter = "top",
                  class    = "table-striped table-hover table-sm")
  })

  observeEvent(input$apply_palette, {
    rv$active_palette <- input$net_palette
    rv$active_dir_pal <- input$edge_dir_pal
    rv$palette_rev    <- FALSE   # reset swap on new palette apply
  })

  observeEvent(input$reverse_palette, {
    rv$palette_rev <- !rv$palette_rev
  })

  # ── Check if network has both positive and negative edges ────────────────────
  has_mixed_signs <- reactive({
    check_neg <- function(df) {
      !is.null(df) && ncol(df) >= 3L && any(as.integer(df[[3L]]) < 0L, na.rm = TRUE)
    }
    check_neg(rv$within_net) || check_neg(rv$btw_net)
  })

  output$edge_sign_toggle_ui <- renderUI({
    req(rv$result)
    if (!has_mixed_signs()) return(NULL)
    div(class = "mt-1",
        checkboxInput("show_edge_sign",
                      "Dashed = negative-correlation edges",
                      value = TRUE))
  })

  # ── Interactive network (visNetwork) — auto-renders after analysis; button freezes to a chosen threshold ─────────
  output$interactive_net <- visNetwork::renderVisNetwork({
    req(rv$result, rv$result$EdgesSelected)
    # Wait for net_group_select_ui to finish rendering so all downstream inputs
    # are ready — prevents 3-4 rapid re-renders from input null-transitions
    req(!is.null(input$net_group_color))
    fe <- if (isTRUE(input$apply_threshold > 0)) filtered_edges_applied() else filtered_edges()
    group_sel   <- input$net_group_color  # reactive dep: re-renders on group change
    if (is.null(fe) || nrow(fe) == 0) return(NULL)

    applied_lambda <- isolate(input$thres_slider)
    if (is.null(applied_lambda)) applied_lambda <- 0

    # ── Nodes ──────────────────────────────────────────────────────────────────
    all_nodes  <- unique(c(fe$NodeA, fe$NodeB))
    data_type  <- ifelse(grepl("datY", all_nodes, fixed = TRUE), "Data Y", "Data X")
    pal_map <- list(
      classic = c(X = "#1a6fafB3", Y = "#e67e22B3"),
      purple  = c(X = "#2471A3B3", Y = "#8E44ADB3"),
      cool    = c(X = "#16a085B3", Y = "#c0392bB3"),
      mono    = c(X = "#555555B3", Y = "#aaaaaa80")
    )
    dir_pal_map <- list(
      rb  = c(up = "#e74c3c", down = "#3498db"),
      br  = c(up = "#3498db", down = "#e74c3c"),
      gr  = c(up = "#27ae60", down = "#e74c3c"),
      bg  = c(up = "#2980b9", down = "#f39c12"),
      tg  = c(up = "#00897b", down = "#1565c0")
    )
    pal_key    <- if (!is.null(rv$active_palette) && rv$active_palette %in% names(pal_map))
                   rv$active_palette else "classic"
    pal        <- pal_map[[pal_key]]
    if (isTRUE(rv$palette_rev)) pal <- c(X = pal[["Y"]], Y = pal[["X"]])
    node_color <- ifelse(data_type == "Data Y", pal[["Y"]], pal[["X"]])
    node_label <- sub("\\.\\.\\.(datX|datY)$", "", all_nodes)
    deg        <- table(c(fe$NodeA, fe$NodeB))
    node_deg   <- as.integer(deg[all_nodes])
    node_deg[is.na(node_deg)] <- 1L
    deg_range  <- range(node_deg)
    node_size  <- if (deg_range[1L] == deg_range[2L]) rep(80L, length(node_deg)) else
      40L + as.integer(120L * (node_deg - deg_range[1L]) / (deg_range[2L] - deg_range[1L]))
    node_font_size <- if (deg_range[1L] == deg_range[2L]) rep(40L, length(node_deg)) else
      28L + as.integer(52L * (node_deg - deg_range[1L]) / (deg_range[2L] - deg_range[1L]))

    # Layout: reuse cached layout if the node set hasn't changed (preserves
    # positions when only edge colors/style change), else recompute.
    layout_key <- paste(sort(all_nodes), collapse = "\x1f")
    cached_lo    <- isolate(rv$net_layout)
    cached_key   <- isolate(rv$net_layout_key)
    cached_nodes <- isolate(rv$net_layout_nodes)
    if (!is.null(cached_lo) && identical(cached_key, layout_key) && !is.null(cached_nodes)) {
      node_idx <- match(all_nodes, cached_nodes)
      lo <- cached_lo[node_idx, , drop = FALSE]
    } else {
      g_tmp <- igraph::graph_from_data_frame(
        data.frame(from = fe$NodeA, to = fe$NodeB, stringsAsFactors = FALSE),
        directed = FALSE,
        vertices = data.frame(name = all_nodes, stringsAsFactors = FALSE)
      )
      lo <- igraph::layout_with_fr(g_tmp, niter = 300) * 200L
      isolate({
        rv$net_layout       <- lo
        rv$net_layout_key   <- layout_key
        rv$net_layout_nodes <- all_nodes
      })
    }
    nodes_df <- data.frame(
      id        = all_nodes, label = node_label, group = data_type,
      color     = node_color, size  = node_size,
      font.size = node_font_size,
      x         = lo[, 1],
      y         = lo[, 2],
      title     = paste0(node_label, " (", data_type, ") | degree: ", node_deg),
      stringsAsFactors = FALSE
    )

    # ── Sign lookup (used for both dash pattern and magnitude edge coloring) ─────
    sign_lookup <- c()
    for (net_df in list(rv$within_net, rv$btw_net)) {
      if (!is.null(net_df) && ncol(net_df) >= 3L) {
        a <- as.character(net_df[[1L]]); b <- as.character(net_df[[2L]])
        keys <- paste(pmin(a, b), pmax(a, b), sep = "___")
        vals <- setNames(as.integer(net_df[[3L]]), keys)
        sign_lookup <- c(sign_lookup, vals)
      }
    }
    na_clean <- sub("\\.\\.\\.(datX|datY)$", "", fe$NodeA)
    nb_clean <- sub("\\.\\.\\.(datX|datY)$", "", fe$NodeB)
    ekeys    <- paste(pmin(na_clean, nb_clean), pmax(na_clean, nb_clean), sep = "___")
    esigns   <- sign_lookup[ekeys]
    esigns[is.na(esigns)] <- 1L

    # ── Edges: magnitude mode vs group-direction mode ──────────────────────────
    fe_cn    <- names(fe); if (is.null(fe_cn)) fe_cn <- colnames(fe)
    abs_cols <- grep("_absdiknew$", fe_cn, value = TRUE, ignore.case = TRUE)
    use_group <- !is.null(group_sel) && group_sel != "magnitude" && nchar(group_sel) > 0

    if (use_group) {
      dir_col <- paste0(group_sel, "_dir")
      mag_col <- paste0(group_sel, "_absdiknew")
      has_cols <- dir_col %in% fe_cn && mag_col %in% fe_cn
    } else {
      has_cols <- FALSE
    }

    dir_pal_key <- if (!is.null(rv$active_dir_pal) && rv$active_dir_pal %in% names(dir_pal_map))
                    rv$active_dir_pal else "rb"
    dir_pal     <- dir_pal_map[[dir_pal_key]]

    if (use_group && has_cols) {
      edge_colors <- ifelse(fe[[dir_col]] == "up",   dir_pal[["up"]],
                     ifelse(fe[[dir_col]] == "down",  dir_pal[["down"]], "#aaaaaa"))
      mag_vals  <- as.numeric(fe[[mag_col]])
      mag_vals[is.na(mag_vals)] <- 0
      mag_min  <- min(mag_vals, na.rm = TRUE)
      mag_max  <- max(mag_vals, na.rm = TRUE)
      mag_norm <- if (mag_max > mag_min) (mag_vals - mag_min) / (mag_max - mag_min) else rep(0.5, length(mag_vals))
      edge_title <- sprintf("Edge: %s–%s | %s |dik*| = %.3f | dir: %s",
                            fe$GeneA, fe$GeneB, group_sel, mag_vals, fe[[dir_col]])
      legend_edges <- data.frame(label = c("Up", "Down"), color = c(dir_pal[["up"]], dir_pal[["down"]]),
                                 stringsAsFactors = FALSE)
    } else {
      # Edge magnitude: weight from network partial correlations (abs_cols), default 1 if absent
      if (length(abs_cols) > 0) {
        avg_mag <- rowMeans(fe[, abs_cols, drop = FALSE], na.rm = TRUE)
      } else {
        avg_mag <- rep(1, nrow(fe))
      }
      mag_min  <- min(avg_mag, na.rm = TRUE)
      mag_max  <- max(avg_mag, na.rm = TRUE)
      mag_norm <- if (mag_max > mag_min) (avg_mag - mag_min) / (mag_max - mag_min) else rep(0.5, length(avg_mag))
      edge_colors <- ifelse(esigns >= 0L, dir_pal[["up"]], dir_pal[["down"]])
      edge_title  <- sprintf("Edge: %s–%s | avg |dik*| = %.3f | sign: %s",
                             fe$GeneA, fe$GeneB, avg_mag,
                             ifelse(esigns >= 0L, "positive", "negative"))
      legend_edges <- data.frame(label = c("Positive", "Negative"),
                                 color = c(dir_pal[["up"]], dir_pal[["down"]]),
                                 stringsAsFactors = FALSE)
    }

    edge_dashes <- esigns < 0L
    # Only show dashes when network has mixed signs AND the toggle is on
    use_dash <- any(edge_dashes) && any(!edge_dashes) && !isFALSE(input$show_edge_sign)
    dash_val <- lapply(edge_dashes, function(d) if (use_dash && isTRUE(d)) c(34L, 20L) else FALSE)
    edges_df <- data.frame(
      from  = fe$NodeA, to = fe$NodeB,
      color = edge_colors, title = edge_title,
      width = 3 + 30 * mag_norm,
      stringsAsFactors = FALSE
    )
    edges_df$dashes <- dash_val

    net <- visNetwork::visNetwork(
      nodes_df, edges_df,
      main = paste0("Subnetwork at λ = ", round(applied_lambda, 2),
                    " (", nrow(fe), " edges)")
    ) %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = list(enabled = TRUE, useLabels = TRUE,
                                main = "Search / select gene(s)"),
        selectedBy = "group"
      ) %>%
      visNetwork::visLegend(
        addNodes = data.frame(
          label = c("Data X", "Data Y"),
          color = c(substring(pal[["X"]], 1, 7), substring(pal[["Y"]], 1, 7)),
          shape = "dot", size = 20, stringsAsFactors = FALSE),
        addEdges = if (!is.null(legend_edges))
          data.frame(label = legend_edges$label, color = legend_edges$color,
                     stringsAsFactors = FALSE) else NULL,
        useGroups = FALSE, position = "right", width = 0.15
      ) %>%
      visNetwork::visPhysics(enabled = FALSE) %>%
      visNetwork::visNodes(
        font    = list(strokeWidth = 3, strokeColor = "#ffffff", color = "#333333"),
        scaling = list(label = list(enabled = FALSE, drawThreshold = 1))
      ) %>%
      visNetwork::visEdges(smooth = if (isTRUE((input$edge_smooth_type %||% "continuous") == "none"))
                                      list(enabled = FALSE)
                                    else
                                      list(type = input$edge_smooth_type %||% "continuous")) %>%
      visNetwork::visInteraction(navigationButtons = TRUE) %>%
      visNetwork::visEvents(
        dragEnd = htmlwidgets::JS(
          "function(params) {
             if (params.nodes && params.nodes.length > 0 && window.Shiny) {
               Shiny.setInputValue('interactive_net_positions',
                 this.getPositions(), {priority: 'event'});
             }
           }"
        )
      )
    net
  })

  # Capture user-moved node positions → update layout cache
  observeEvent(input$interactive_net_positions, {
    pos <- input$interactive_net_positions
    if (is.null(pos) || length(pos) == 0L || is.null(rv$net_layout_nodes)) return()
    nodes  <- rv$net_layout_nodes
    lo_new <- matrix(0, nrow = length(nodes), ncol = 2L)
    for (i in seq_along(nodes)) {
      nd_pos <- pos[[as.character(nodes[i])]]
      if (!is.null(nd_pos)) {
        lo_new[i, 1L] <- if (!is.null(nd_pos$x)) as.numeric(nd_pos$x) else 0
        lo_new[i, 2L] <- if (!is.null(nd_pos$y)) as.numeric(nd_pos$y) else 0
      } else if (!is.null(rv$net_layout) && i <= nrow(rv$net_layout)) {
        lo_new[i, ] <- rv$net_layout[i, ]
      }
    }
    rv$net_layout <- lo_new
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # ── Network hub table ──────────────────────────────────────────────────────
  output$net_hub_tbl <- DT::renderDT({
    fe <- if (isTRUE(input$apply_threshold > 0)) filtered_edges_applied() else
          if (!is.null(rv$result) && !is.null(rv$result$EdgesSelected)) rv$result$EdgesSelected else NULL
    if (is.null(fe)) return(DT::datatable(data.frame()))
    deg <- sort(table(c(fe$NodeA, fe$NodeB)), decreasing = TRUE)
    top <- head(deg, 10)
    data.frame(
      Node   = sub("\\.\\.\\.(datX|datY)$", "", names(top)),
      Type   = ifelse(grepl("datY", names(top), fixed = TRUE), "Data Y", "Data X"),
      Degree = as.integer(top),
      stringsAsFactors = FALSE
    )
  }, options = list(dom = "t", pageLength = 10, ordering = FALSE),
     rownames = FALSE, selection = "none")

  # ── Network top-edge table ─────────────────────────────────────────────────
  output$net_edge_tbl_hdr <- renderUI({
    grp <- input$net_group_color %||% "magnitude"
    lbl <- if (grp == "magnitude") "Top 10 Edges (avg dik*)" else
             paste0("Top 10 Edges – ", grp, " dik*")
    h5(class = "mb-2 mt-2", lbl)
  })

  output$net_edge_tbl <- DT::renderDT({
    fe <- if (isTRUE(input$apply_threshold > 0)) filtered_edges_applied() else
          if (!is.null(rv$result) && !is.null(rv$result$EdgesSelected)) rv$result$EdgesSelected else NULL
    if (is.null(fe)) return(DT::datatable(data.frame()))
    fe_cn  <- names(fe); if (is.null(fe_cn)) fe_cn <- colnames(fe)
    grp    <- input$net_group_color %||% "magnitude"
    if (grp == "magnitude") {
      abs_cols <- grep("_absdiknew$", fe_cn, value = TRUE, ignore.case = TRUE)
      if (length(abs_cols) == 0) return(DT::datatable(data.frame()))
      # sort by absolute magnitude; display signed diknew
      sort_vals  <- rowMeans(fe[, abs_cols, drop = FALSE], na.rm = TRUE)
      dik_cols   <- sub("(?i)_absdiknew$", "_diknew", abs_cols, perl = TRUE)
      dik_cols   <- dik_cols[dik_cols %in% fe_cn]
      score_vals <- if (length(dik_cols) > 0)
        rowMeans(fe[, dik_cols, drop = FALSE], na.rm = TRUE) else sort_vals
      col_lbl    <- "avg dik*"
    } else {
      abs_col  <- paste0(grp, "_absdiknew")
      dik_col  <- paste0(grp, "_diknew")
      if (!abs_col %in% fe_cn) return(DT::datatable(data.frame()))
      sort_vals  <- as.numeric(fe[[abs_col]])
      score_vals <- if (dik_col %in% fe_cn) as.numeric(fe[[dik_col]]) else sort_vals
      col_lbl    <- paste0(grp, " dik*")
    }
    top_idx <- order(sort_vals, decreasing = TRUE)[seq_len(min(10L, nrow(fe)))]
    top_fe  <- fe[top_idx, ]
    out <- data.frame(
      GeneA  = top_fe$GeneA,
      GeneB  = top_fe$GeneB,
      Score  = round(score_vals[top_idx], 3),
      stringsAsFactors = FALSE
    )
    if (grp != "magnitude" && paste0(grp, "_dir") %in% fe_cn)
      out$Direction <- top_fe[[paste0(grp, "_dir")]]
    names(out)[3] <- col_lbl
    out
  }, options = list(dom = "t", pageLength = 10, ordering = FALSE),
     rownames = FALSE, selection = "none")

  # ── ZIP download ───────────────────────────────────────────────────────────
  output$dl_all_zip <- downloadHandler(
    filename = function() {
      paste0("iOmicsPASS_results_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(zipfile) {
      req(rv$result)
      zip_dir  <- file.path(tempdir(), "iOmicsPASS_zip_export")
      unlink(zip_dir, recursive = TRUE)
      dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)

      write_if <- function(obj, fname) {
        if (is.null(obj)) return(invisible(NULL))
        if (is.data.frame(obj) || is.matrix(obj))
          write.table(obj, file.path(zip_dir, fname),
                      sep = "\t", quote = FALSE, row.names = FALSE)
        else if (is.character(obj))
          writeLines(obj, file.path(zip_dir, fname))
      }

      write_if(rv$result$EdgesSelected,       "EdgesSelected.txt")
      write_if(rv$result$SampleProbabilities, "SampleProbabilities.txt")
      write_if(rv$result$CVerrors,            "CVerrors.txt")
      write_if(rv$result$AttributesTable,     "AttributesTable.txt")
      write_if(rv$result$BGlist,              "BGlist.txt")
      write_if(rv$result$ZtransX,             "ZtransX.txt")
      write_if(rv$result$ZtransY,             "ZtransY.txt")
      write_if(rv$result$ZtransZ,             "ZtransZ.txt")

      # Gene-level pathway enrichment (if computed)
      gpd <- tryCatch(isolate(gene_pathway_data()), error = function(e) NULL)
      if (!is.null(gpd) && nrow(gpd) > 0)
        write.table(gpd, file.path(zip_dir, "GenePathwayEnrichment.txt"),
                    sep = "\t", row.names = FALSE, quote = FALSE)

      res <- rv$result
      summary_lines <- c(
        paste0("iOmicsPASS+ Analysis Summary — ", Sys.time()),
        paste0("Optimal threshold (minThres): ", res$minThres),
        paste0("Selected edges: ",
               if (!is.null(res$EdgesSelected)) nrow(res$EdgesSelected) else "NA")
      )
      writeLines(summary_lines, file.path(zip_dir, "summary.txt"))

      # Pull any extra files written to tmpDir by the Rcpp core.
      # Skip Ztransform_data*.txt (already included as ZtransX/Y/Z.txt above)
      # and SampleClass_Probabilities.txt (already included as SampleProbabilities.txt above).
      if (!is.null(res$tmpDir) && dir.exists(res$tmpDir)) {
        extra <- list.files(res$tmpDir, full.names = TRUE, recursive = FALSE)
        for (f in extra) {
          if (grepl("^Ztransform_data", basename(f))) next
          if (basename(f) == "SampleClass_Probabilities.txt") next
          dest <- file.path(zip_dir, basename(f))
          if (!file.exists(dest)) file.copy(f, dest)
        }
      }

      # ── Plots ──────────────────────────────────────────────────────────────
      # Network plot (captured from last Save PNG/JPEG in Network Summary tab)
      if (!is.null(rv$net_export_path) && file.exists(rv$net_export_path))
        file.copy(rv$net_export_path, file.path(zip_dir, "Network_plot.png"))

      # CV plot
      tryCatch({
        cv_png <- file.path(zip_dir, "CVplot.png")
        png(cv_png, width = 1500, height = 1000, res = 200)
        .render_cv_plot()
        dev.off()
      }, error = function(e) NULL)

      # Network enrichment bubble plot
      tryCatch({
        ned  <- isolate(pathway_data())
        if (!is.null(ned) && nrow(ned) > 0) {
          cls  <- isolate(input$pathway_class) %||% ""
          pcut <- isolate(input$pathway_pval)  %||% 0.05
          dims <- .bubble_dims(ned)
          ne_png <- file.path(zip_dir, "NetworkEnrichment_plot.png")
          png(ne_png, width = dims$width, height = dims$height, res = 200)
          .draw_net_enrich_plot(ned, cls, pcut, left_mar = dims$left_mar)
          dev.off()
        }
      }, error = function(e) NULL)

      # Gene pathway enrichment bubble plot
      tryCatch({
        ged  <- isolate(gene_pathway_data())
        if (!is.null(ged) && nrow(ged) > 0) {
          grp  <- isolate(input$gene_enrich_group) %||% ""
          pcut <- isolate(input$gene_enrich_pval)  %||% 0.05
          dims <- .bubble_dims(ged)
          ge_png <- file.path(zip_dir, "PathwayEnrichment_plot.png")
          png(ge_png, width = dims$width, height = dims$height, res = 200)
          .draw_gene_enrich_plot(ged, grp, pcut, left_mar = dims$left_mar)
          dev.off()
        }
      }, error = function(e) NULL)

      all_files <- list.files(zip_dir, full.names = TRUE)
      owd <- getwd(); setwd(zip_dir); on.exit(setwd(owd), add = TRUE)
      zip(zipfile, basename(all_files))
    }
  )

  # ── Prediction tab ─────────────────────────────────────────────────────────
  output$predict_model_status <- renderUI({
    if (is.null(rv$result) || is.null(rv$result$EdgesSelected)) {
      div(class = "alert alert-warning",
          icon("exclamation-triangle"),
          " No analysis results available. ",
          "Please complete the iOmicsPASS Analysis first.")
    } else {
      div(class = "alert alert-success",
          icon("check-circle"),
          sprintf(" Model ready: %d selected edges available for prediction.",
                  nrow(rv$result$EdgesSelected)))
    }
  })

  observeEvent(input$run_predict, {
    req(rv$result, rv$result$EdgesSelected)
    req(input$file_newdata)

    withProgress(message = "Running prediction on new samples...", value = 0.5, {
      tryCatch({
        new_data <- read_upload(input$file_newdata$datapath)
        if (is.null(new_data)) stop("Could not read new data file.")

        pred <- Predict.iOmicsPASS(
          file        = rv$result$EdgesSelected,
          newData     = new_data,
          standardize = isTRUE(input$predict_standardize),
          prop        = input$predict_prop
        )

        rv$pred_result <- pred
        showNotification(
          sprintf("Prediction complete: %d samples classified.", nrow(pred)),
          type = "message", duration = 3L
        )
      }, error = function(e) {
        showNotification(paste("Prediction error:", conditionMessage(e)),
                         type = "error", duration = 8L)
      })
    })
  })

  output$predict_results_ui <- renderUI({
    if (is.null(rv$pred_result)) {
      return(div(class = "alert alert-info",
                 icon("info-circle"),
                 " Upload new data and click Run Prediction."))
    }
    tagList(
      div(class = "alert alert-success",
          icon("check-circle"),
          sprintf(" %d sample(s) classified.", nrow(rv$pred_result))),
      downloadButton("dl_predictions", "Download predictions (.txt)",
                     class = "btn btn-sm btn-outline-primary mb-2")
    )
  })

  output$tbl_predict <- DT::renderDataTable({
    if (is.null(rv$pred_result)) return(DT::datatable(data.frame()))
    DT::datatable(rv$pred_result,
                  options   = list(pageLength = 10L, scrollX = TRUE),
                  rownames  = FALSE,
                  class     = "table-striped table-hover table-sm")
  })

  output$dl_predictions <- downloadHandler(
    filename = function() paste0("Predictions_", Sys.Date(), ".txt"),
    content  = function(file) {
      write.table(rv$pred_result, file,
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
