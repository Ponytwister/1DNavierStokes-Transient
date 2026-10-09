#include <sqlite3.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// CTest helper: creates synthetic input or verifies the CLI's persisted output.
int main(int argc, char** argv)
{
    try {
        if (argc < 4 || argc > 5) throw std::runtime_error("Expected mode, database, fixture/expected-value, and optional save count");
        const std::string mode = argv[1];
        sqlite3* raw = nullptr;
        const int rc = sqlite3_open_v2(argv[2], &raw,
            mode == "create" ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE : SQLITE_OPEN_READONLY, nullptr);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        if (rc != SQLITE_OK) throw std::runtime_error("Cannot open disposable database");
        if (mode == "create") {
            if (argc != 4) throw std::runtime_error("Create mode expects one fixture path");
            std::ifstream fixture(argv[3]);
            if (!fixture) throw std::runtime_error("Cannot read fixture");
            const std::string sql(std::istreambuf_iterator<char>(fixture), {});
            if (sqlite3_exec(db.get(), sql.c_str(), nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db.get()));
            const auto migration_path = std::filesystem::path(argv[3]).parent_path() /
                "../../migrations/007_profile_history_natural_keys.sql";
            std::ifstream migration(migration_path);
            if (!migration) throw std::runtime_error("Cannot read profile history migration");
            const std::string migration_sql(std::istreambuf_iterator<char>(migration), {});
            if (sqlite3_exec(db.get(), migration_sql.c_str(), nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db.get()));
            const char* reseed =
                "INSERT INTO raw_profile_entrance_concentrations "
                "SELECT NAME,WT_PERCENT,1,'FITC',0.002,'mg/ml' FROM raw_profile;"
                "INSERT INTO profile_calculations VALUES(99,'unrelated','keq1','all reactions enabled','none',8,4);"
                "INSERT INTO profile_records VALUES(99,'unrelated','99',3,4,5,6);"
                "INSERT INTO profile_record_parameters VALUES(99,'','','keq1',42,'','global');"
                "INSERT INTO profile_record_samples(RUN_AT_UTC_US,EXPERIMENT_NAME,WT_PERCENT,X,\"Numeric\") "
                "VALUES(99,'unrelated','99',0,42);";
            if (sqlite3_exec(db.get(), reseed, nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db.get()));
        } else if (mode == "check") {
            if (argc != 5) throw std::runtime_error("Check mode expects fitted value and profile save count");
            const int profile_saves = std::stoi(argv[4]);
            if (profile_saves < 1 || profile_saves > 2) throw std::runtime_error("Unexpected profile save count");
            const std::string sql =
                "SELECT (SELECT count(*) FROM profile_record_samples WHERE RUN_AT_UTC_US<>99)=" + std::to_string(36 * profile_saves) +
                " AND (SELECT count(*) FROM profile_record_samples WHERE RUN_AT_UTC_US<>99 AND abs(\"Numeric\"-0.002)<=1e-12)=" + std::to_string(36 * profile_saves) +
                " AND (SELECT count(*) FROM profile_record_parameters WHERE RUN_AT_UTC_US<>99)=" + std::to_string(profile_saves) +
                " AND (SELECT count(*) FROM profile_calculations)=" + std::to_string(1 + profile_saves) +
                " AND (SELECT count(*) FROM profile_records)=" + std::to_string(1 + 4 * profile_saves) +
                " AND (SELECT \"Numeric\" FROM profile_record_samples WHERE RUN_AT_UTC_US=99)=42"
                " AND (SELECT [INITIAL VALUE] FROM alglib_input)=" + std::string(argv[3]);
            sqlite3_stmt* stmt = nullptr;
            const int rc = sqlite3_prepare_v2(db.get(), sql.c_str(), -1, &stmt, nullptr);
            std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(stmt, sqlite3_finalize);
            if (rc != SQLITE_OK || sqlite3_step(stmt) != SQLITE_ROW || sqlite3_column_int(stmt, 0) != 1)
                throw std::runtime_error("CLI persistence did not match the independent fixture");
        } else throw std::runtime_error("Unknown mode");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
