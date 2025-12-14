using Dates
using Compose

# Paths
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT_PATH = "./results-02/"

# Scenario defaults
const DEFAULT_PATIENT_TYPE = :icu
const CAPACITY_LEVELS = ["baseline", "rampup", "surge", "nextup", "max"]
const CAPACITY_NAMES = ["Baseline", "Level 1", "Level 2", "Level 3", "Max"]

# Plot styling
const DEFAULT_FONT = "Helvetica"
const FONTSTYLES = (
    key_title_font = DEFAULT_FONT,
    key_label_font = DEFAULT_FONT,
    minor_label_font = DEFAULT_FONT,
    major_label_font = DEFAULT_FONT,
    point_label_font = DEFAULT_FONT,
    key_title_font_size = 14px,
    key_label_font_size = 13px,
    minor_label_font_size = 14px,
    major_label_font_size = 16px,
    point_label_font_size = 12px,
)

halfmonth(d::Date) = dayofmonth(d) < 15 ? firstdayofmonth(d) : firstdayofmonth(d) + Day(14)

function date_ticks_for_range(start_date::Date, end_date::Date)
    ds = collect(start_date:Day(1):end_date)
    ds = unique(halfmonth.(ds))
    return DateTime.(ds)
end
