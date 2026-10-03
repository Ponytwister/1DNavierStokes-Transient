#include <tsensor.h>
#include <memory>
#include <array>

void initialize_channel_dimensions(parameters_t& p, sqlite3* db) {
    if (!p.experiments.empty()) throw std::logic_error("Channel geometry must be initialized before model links.");
    std::vector<experiment_run_struct> runs;
    runs.reserve(p.experiment_runs.size());
    for (const auto& selected : p.experiment_runs) {
        check_cancellation(p);
        sqlite3_stmt* raw = nullptr;
        const int prepared = sqlite3_prepare_v2(db, "SELECT * FROM experiments WHERE NAME=?", -1, &raw, nullptr);
        std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(raw, sqlite3_finalize);
        if (prepared != SQLITE_OK) throw std::runtime_error(sqlite3_errmsg(db));
        if (sqlite3_bind_text(raw, 1, selected.name.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK)
            throw std::runtime_error(sqlite3_errmsg(db));
        const int first = sqlite3_step(raw);
        if (first != SQLITE_ROW) throw std::runtime_error("Cannot load channel geometry for experiment: " + selected.name);
        std::array<int, 3> columns{-1, -1, -1};
        const char* names[] = {"CHANNEL_WIDTH", "CHANNEL_HEIGHT", "CHANNEL_LENGTH"};
        for (int col = 0; col < sqlite3_column_count(raw); ++col)
            for (int i = 0; i < 3; ++i)
                if (sqlite3_stricmp(sqlite3_column_name(raw, col), names[i]) == 0) columns[i] = col;
        const int present = std::count_if(columns.begin(), columns.end(), [](int col) { return col >= 0; });
        if (present != 0 && present != 3) throw std::invalid_argument("Incomplete channel schema; apply the channel dimensions migration.");
        double values[] = {p.W, p.H, p.L};
        if (present == 3) {
            for (int i = 0; i < 3; ++i) {
                const int type = sqlite3_column_type(raw, columns[i]);
                if (type != SQLITE_FLOAT && type != SQLITE_INTEGER)
                    throw std::invalid_argument(selected.name + ": " + names[i] + " must be numeric and non-NULL.");
                values[i] = sqlite3_column_double(raw, columns[i]);
            }
        }
        const channel_dimensions dimensions(values[0], values[1], values[2]);
        const int next = sqlite3_step(raw);
        if (next != SQLITE_DONE) throw std::runtime_error("Ambiguous or unreadable experiment: " + selected.name);
        runs.emplace_back(dimensions);
        runs.back().name = selected.name;
    }
    p.experiment_runs.swap(runs);
}
