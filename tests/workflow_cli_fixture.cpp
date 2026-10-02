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
        if (argc != 4) throw std::runtime_error("Expected mode, database, fixture/expected-value");
        const std::string mode = argv[1];
        sqlite3* raw = nullptr;
        const int rc = sqlite3_open_v2(argv[2], &raw,
            mode == "create" ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE : SQLITE_OPEN_READONLY, nullptr);
        std::unique_ptr<sqlite3, decltype(&sqlite3_close)> db(raw, sqlite3_close);
        if (rc != SQLITE_OK) throw std::runtime_error("Cannot open disposable database");
        if (mode == "create") {
            std::ifstream fixture(argv[3]);
            if (!fixture) throw std::runtime_error("Cannot read fixture");
            const std::string sql(std::istreambuf_iterator<char>(fixture), {});
            if (sqlite3_exec(db.get(), sql.c_str(), nullptr, nullptr, nullptr) != SQLITE_OK)
                throw std::runtime_error(sqlite3_errmsg(db.get()));
        } else if (mode == "check") {
            const std::string sql =
                "SELECT (SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99)=36"
                " AND (SELECT count(*) FROM model_profile WHERE SOLUTION_ID<>99 AND abs(Numeric-0.002)<=1e-12)=36"
                " AND (SELECT count(*) FROM parameter_solutions WHERE SOLVE_SETTING_ID<>99)=1"
                " AND (SELECT count(*) FROM solutions)=5"
                " AND (SELECT Numeric FROM model_profile WHERE SOLUTION_ID=99)=42"
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
