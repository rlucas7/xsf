#include "../testing_utils.h"
#include <cmath>
#include <limits>
#include <mdspan>
#include <vector>
#include <xsf/legendre.h>

TEST_CASE("lqmn endpoints return sentinel values", "[lqmn][xsf_tests]") {
    constexpr int m = 3;
    constexpr int n = 4;

    const double x = GENERATE(-1.0, 1.0);

    std::vector<double> qm_buf((m + 1) * (n + 1));
    std::vector<double> qd_buf((m + 1) * (n + 1));

    std::mdspan<double, std::dextents<std::size_t, 2>> qm(qm_buf.data(), m + 1, n + 1);
    std::mdspan<double, std::dextents<std::size_t, 2>> qd(qd_buf.data(), m + 1, n + 1);

    xsf::lqmn(x, qm, qd);

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            CAPTURE(x, i, j, qm(i, j), qd(i, j));
            REQUIRE(qm(i, j) == 1e300);
            REQUIRE(qd(i, j) == 1e300);
        }
    }
}

TEST_CASE("lqmn stays finite next to singular endpoints", "[lqmn][xsf_tests]") {
    constexpr int m = 3;
    constexpr int n = 4;

    const double x = GENERATE(
        std::nextafter(-1.0, -2.0),
        std::nextafter(-1.0, 0.0),
        std::nextafter(1.0, 0.0),
        std::nextafter(1.0, 2.0)
    );

    std::vector<double> qm_buf((m + 1) * (n + 1));
    std::vector<double> qd_buf((m + 1) * (n + 1));

    std::mdspan<double, std::dextents<std::size_t, 2>> qm(qm_buf.data(), m + 1, n + 1);
    std::mdspan<double, std::dextents<std::size_t, 2>> qd(qd_buf.data(), m + 1, n + 1);

    xsf::lqmn(x, qm, qd);

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            CAPTURE(x, i, j, qm(i, j), qd(i, j));
            REQUIRE(std::isfinite(qm(i, j)));
            REQUIRE(std::isfinite(qd(i, j)));
        }
    }
}
