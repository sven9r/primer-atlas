# A refresh verifies a complete candidate before replacing one session's
# reference. Failed or incomplete latest pointers never select a fallback.
coi_reference_release_identity <- function(reference_or_state) {
  state <- if (!is.null(reference_or_state$state)) reference_or_state$state else reference_or_state
  manifest <- state$manifest
  hashes <- function(artifacts) {
    keys <- sort(names(artifacts))
    setNames(lapply(keys, function(key) artifacts[[key]]$sha256), keys)
  }
  retained <- manifest$coi_retained_state
  value <- list(marker_id = manifest$marker_id, release_version = manifest$release_version,
                artifacts = hashes(manifest$artifacts))
  if (!is.null(retained)) {
    value$retained <- retained[c("schema_version", "retained_records", "target_groups", "group_counts", "policy")]
    value$retained$artifacts <- hashes(retained$artifacts)
  }
  digest::digest(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null"), algo = "sha256", serialize = FALSE)
}

coi_reference_check_refresh <- function(previous, candidate, root = getwd()) {
  if (!identical(candidate$mode, "full_release") || !identical(candidate$alignment_mode, "full")) {
    stop("Latest COI release does not supply a complete verified full reference")
  }
  if (identical(previous$alignment_mode, "full") && nrow(previous$records)) {
    if (!setequal(names(previous$counts), names(candidate$counts))) stop("New COI release changes the retained target-group scope")
    before <- previous$records
    after <- candidate$records
    keys <- paste(before$target_group, before$accession)
    matched <- match(keys, paste(after$target_group, after$accession))
    if (anyNA(matched)) stop("New COI release loses previously retained accessions")
    for (column in c("sequence_sha256", "header_sha256")) {
      if (!identical(as.character(before[[column]]), as.character(after[[column]][matched]))) {
        stop("New COI release changes a previously retained ", column)
      }
    }
  }
  # Verify the partition bytes before activation without parsing the large
  # position tables. Their rows remain lazy in memory and are independently
  # checked against the active cohort when a detail view requests them.
  artifacts <- candidate$state$manifest$artifacts
  pairs <- setdiff(unique(candidate$overview$pair_id), "BEEPRIME")
  for (pair in pairs) {
    prefix <- paste0("full_order_", tolower(pair))
    positions <- intersect(c(paste0(prefix, "_positions"), paste0(prefix, "_primer_positions")), names(artifacts))
    if (!prefix %in% names(artifacts) || !length(positions)) stop("New COI release lacks a complete partition for ", pair)
    for (key in c(prefix, positions[1])) {
      artifact <- artifacts[[key]]
      if (is.null(artifact$sha256) || !grepl("^[a-fA-F0-9]{64}$", artifact$sha256) ||
          is.null(artifact$rows) || !is.finite(artifact$rows) || artifact$rows < candidate$n_references) {
        stop("New COI release has an invalid partition descriptor for ", pair)
      }
      coi_reference_artifact_path(artifact, root)
    }
    if (artifacts[[prefix]]$rows != candidate$n_references) stop("New COI sequence partition has a different denominator for ", pair)
  }
  invisible(TRUE)
}

coi_reference_refresh_coordinator <- function(initial, root = getwd(),
    reader = function() atlas_read_manifest("COI", root),
    loader = function(state) .coi_reference_full_release(state, root, 100L, 10000L),
    now = function() format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  current <- initial
  # Caches are mutable; they must never be shared between app sessions.
  current$cache <- new.env(parent = emptyenv())
  current_identity <- coi_reference_release_identity(current)
  status <- list(phase = "ready", checked_at = NULL, message = "Update checks run every 10 minutes.")
  check <- function() {
    candidate <- tryCatch({
      state <- reader()
      if (is.null(state) || isTRUE(state$stale)) stop("The latest COI pointer could not be read from its source")
      hosted <- !is.null(state$source) && grepl("^https?://", state$source)
      source_label <- if (hosted) "hosted" else "bundled"
      candidate_identity <- coi_reference_release_identity(state)
      if (identical(state$manifest$release_version, current$version) && identical(candidate_identity, current_identity)) {
        list(changed = FALSE, reference = current,
          status = list(phase = "current", checked_at = now(), message = if (hosted) {
            "This session has the latest verified hosted reference."
          } else "This session matches the bundled reference. Hosted update checks are not configured."))
      } else {
        if (identical(state$manifest$release_version, current$version)) stop("An immutable COI release changed its content identity")
        if (is.null(state$manifest$coi_retained_state)) stop("The ", source_label,
          " COI reference has no complete growing-reference state.",
          if (!hosted) " Hosted update checks are not configured." else "")
        verified <- loader(state)
        coi_reference_check_refresh(current, verified, root)
        verified$cache <- new.env(parent = emptyenv())
        list(reference = verified, identity = candidate_identity, source_label = source_label)
      }
    }, error = base::identity)
    if (inherits(candidate, "error")) {
      status <<- list(phase = "retained", checked_at = now(),
                     message = paste0("Current reference retained. ", conditionMessage(candidate)))
      return(list(changed = FALSE, reference = current, status = status))
    }
    if (!is.null(candidate$changed)) {
      status <<- candidate$status
      return(candidate)
    }
    previous_version <- current$version
    current <<- candidate$reference
    current_identity <<- candidate$identity
    status <<- list(phase = "updated", checked_at = now(),
                   message = paste0("Verified ", candidate$source_label, " reference updated from ", previous_version, "."))
    list(changed = TRUE, reference = current, status = status)
  }
  list(current = function() current, check = check, status = function() status)
}

coiReferenceRefreshUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(class = "coi-refresh-bar",
    shiny::uiOutput(ns("status")),
    shiny::actionButton(ns("check_latest"), "Check for COI updates", class = "btn-outline-secondary btn-sm"))
}

coiReferenceRefreshServer <- function(id, initial, root = getwd(), interval_ms = 600000L,
    reader = function() atlas_read_manifest("COI", root),
    loader = function(state) .coi_reference_full_release(state, root, 100L, 10000L)) {
  shiny::moduleServer(id, function(input, output, session) {
    coordinator <- coi_reference_refresh_coordinator(initial, root, reader, loader)
    reference <- shiny::reactiveVal(coordinator$current())
    status <- shiny::reactiveVal(coordinator$status())
    check <- function() {
      result <- coordinator$check()
      if (isTRUE(result$changed)) {
        reference(result$reference)
      }
      status(result$status)
      invisible(result)
    }
    shiny::observe({
      shiny::invalidateLater(interval_ms, session)
      shiny::isolate(check())
    })
    shiny::observeEvent(input$check_latest, check(), ignoreInit = TRUE)
    output$status <- shiny::renderUI({
      current <- reference()
      latest <- status()
      shiny::div(
        shiny::strong("COI reference: ", current$version),
        shiny::span(" · ", format(current$n_references, big.mark = ","),
          if (identical(current$alignment_mode, "full")) " retained references" else " published score references · preview alignment only"),
        shiny::div(class = if (identical(latest$phase, "retained")) "tiny text-warning-emphasis" else "tiny",
          latest$message, if (!is.null(latest$checked_at)) paste0(" Checked ", latest$checked_at, ".")))
    })
    list(reference = shiny::reactive(reference()), status = shiny::reactive(status()), check = check,
         identity = shiny::reactive(coi_reference_release_identity(reference())))
  })
}
