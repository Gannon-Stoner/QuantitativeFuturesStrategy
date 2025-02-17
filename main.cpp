#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <curl/curl.h>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// -------------------------------
// Structures and Data Definitions
// -------------------------------

struct Bar {
    std::string timestamp;
    double open;
    double high;
    double low;
    double close;
    long volume;
};

struct Trade {
    std::string entryTime;
    std::string exitTime;
    double entryPrice;
    double exitPrice;
    int direction;  // 1 for long, -1 for short
    double size;
    double profit;
    int holdingPeriod;
};

struct PerformanceMetrics {
    // Return metrics
    double annualizedReturn;
    double annualizedVolatility;
    double sharpeRatio;
    double maxDrawdown;

    // Trade metrics
    int totalTrades;
    int winningTrades;
    int losingTrades;
    double winRate;
    double averageWin;
    double averageLoss;
    double profitFactor;      // Gross profits / Gross losses
    double averageHoldingPeriod;
    double largestWin;
    double largestLoss;
    double maxConsecutiveLosses;
    double calmarRatio;       // Annualized return / Max drawdown
    double sortinoRatio;      // Similar to Sharpe but only considers downside volatility
};

// -------------------------------
// Helper: Curl Write Callback
// -------------------------------
size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    size_t totalSize = size * nmemb;
    std::string* str = static_cast<std::string*>(userp);
    str->append(static_cast<char*>(contents), totalSize);
    return totalSize;
}

// -------------------------------
// DataFetcher Module: Download and Parse Intraday Data
// -------------------------------
std::vector<Bar> fetchHistoricalData(const std::string &symbol,
                                   const std::string &timeframe,
                                   const std::string &fromDate,
                                   const std::string &toDate,
                                   const std::string &apikey) {
    std::vector<Bar> bars;

    // Build URL; for example:
    // https://financialmodelingprep.com/api/v3/historical-chart/5min/ZOUSX?from=2024-02-10&to=2024-03-10&apikey=YOUR_API_KEY_HERE
    std::string url = "https://financialmodelingprep.com/api/v3/historical-chart/" +
                      timeframe + "/" + symbol +
                      "?from=" + fromDate +
                      "&to=" + toDate +
                      "&apikey=" + apikey;

    // Use libcurl to fetch the data.
    CURL *curl = curl_easy_init();
    if (!curl) {
        std::cerr << "Error initializing curl." << std::endl;
        return bars;
    }
    std::string readBuffer;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        std::cerr << "curl_easy_perform() failed: "
                  << curl_easy_strerror(res) << std::endl;
        curl_easy_cleanup(curl);
        return bars;
    }
    curl_easy_cleanup(curl);

    // Parse the JSON response.
    try {
        auto j = json::parse(readBuffer);
        if (!j.is_array()) {
            std::cerr << "Expected an array of observations for symbol " << symbol << std::endl;
            std::cerr << "Raw JSON: " << j.dump(2) << std::endl;
            return bars;
        }
        for (auto& obs : j) {
            Bar bar;
            bar.timestamp = obs["date"].get<std::string>();
            bar.open = obs["open"].get<double>();
            bar.high = obs["high"].get<double>();
            bar.low  = obs["low"].get<double>();
            bar.close = obs["close"].get<double>();
            bar.volume = obs["volume"].get<long>();
            bars.push_back(bar);
        }
        // Sort bars in chronological order.
        std::sort(bars.begin(), bars.end(), [](const Bar &a, const Bar &b) {
            return a.timestamp < b.timestamp;
        });
    }
    catch (std::exception &e) {
        std::cerr << "JSON parse error for symbol " << symbol << ": " << e.what() << std::endl;
    }

    return bars;
}

// -------------------------------
// Indicator Functions
// -------------------------------

double computeSMA(const std::vector<Bar>& bars, int endIndex, int window) {
    if (endIndex - window + 1 < 0) return 0.0;
    double sum = 0.0;
    for (int i = endIndex - window + 1; i <= endIndex; i++) {
        sum += bars[i].close;
    }
    return sum / window;
}

double computeEMA(const std::vector<Bar>& bars, int endIndex, int window) {
    if (endIndex - window + 1 < 0) return bars[endIndex].close;
    double alpha = 2.0 / (window + 1);
    double ema = bars[endIndex - window + 1].close;
    for (int i = endIndex - window + 2; i <= endIndex; i++) {
        ema = (bars[i].close - ema) * alpha + ema;
    }
    return ema;
}

double computeStdDev(const std::vector<Bar>& bars, int endIndex, int window) {
    if (endIndex - window + 1 < 0) return 0.0;
    double mean = computeSMA(bars, endIndex, window);
    double variance = 0.0;
    for (int i = endIndex - window + 1; i <= endIndex; i++) {
        variance += std::pow(bars[i].close - mean, 2);
    }
    return std::sqrt(variance / window);
}

double zScore(double currentPrice, double mean, double stdDev) {
    if (stdDev == 0) return 0.0;
    return (currentPrice - mean) / stdDev;
}

double computeEMA_ATR(const std::vector<Bar>& bars, int endIndex, int window) {
    if (endIndex - window + 1 < 1) return 0.0;

    double firstTR = std::max({
        bars[endIndex - window + 1].high - bars[endIndex - window + 1].low,
        std::abs(bars[endIndex - window + 1].high - bars[endIndex - window].close),
        std::abs(bars[endIndex - window + 1].low - bars[endIndex - window].close)
    });

    double ema_atr = firstTR;
    double alpha = 2.0 / (window + 1);

    for (int i = endIndex - window + 2; i <= endIndex; i++) {
        double tr = std::max({
            bars[i].high - bars[i].low,
            std::abs(bars[i].high - bars[i-1].close),
            std::abs(bars[i].low - bars[i-1].close)
        });
        ema_atr = (tr - ema_atr) * alpha + ema_atr;
    }

    return ema_atr;
}

int momentumSignal(double currentPrice, double shortEMA, double mediumEMA) {
    if (shortEMA > mediumEMA && currentPrice > shortEMA) return 1;
    if (shortEMA < mediumEMA && currentPrice < shortEMA) return -1;
    return 0;
}

double reactionMultiplier(double currentPrice, double longMA, double stdDev) {
    double z = zScore(currentPrice, longMA, stdDev);
    double x = 1.0 - std::abs(z) / 3.0;  // More tolerant of deviations
    if (x < 0) return 0.0;

    double baseMultiplier = x * std::exp(0.2 * (1 - x * x));

    // Asymmetric adjustment for both directions
    if (currentPrice > longMA || currentPrice < longMA) {
        return baseMultiplier * 0.9;  // 15% reduction for extended moves in either direction
    }

    return baseMultiplier;
}

// -------------------------------
// Backtest Engine
// -------------------------------
struct BacktestResult {
    std::vector<std::string> timestamps;
    std::vector<double> portfolioValues;
    std::vector<double> returns;
    std::vector<Trade> trades;  // New: track all trades
};

BacktestResult runBacktest(const std::vector<Bar>& bars, double initialCapital) {
    BacktestResult result;

    const int shortEmaWindow = 20;
    const int mediumEmaWindow = 50;
    const int longEmaWindow = 100;
    const int atrWindow = 10;
    const double riskUnit = 1000.0;

    int startIndex = longEmaWindow;
    if ((size_t)startIndex >= bars.size()) {
        std::cerr << "Not enough bars for backtesting." << std::endl;
        return result;
    }

    double portfolio = initialCapital;
    bool inTrade = false;
    int tradeDirection = 0;
    double entryPrice = 0.0;
    double positionSize = 0.0;
    int entryIndex = 0;  // Track entry index for holding period

    for (size_t i = startIndex; i < bars.size(); i++) {
        double currentPrice = bars[i].close;
        double shortEMA = computeEMA(bars, i, shortEmaWindow);
        double mediumEMA = computeEMA(bars, i, mediumEmaWindow);
        double longEMA = computeEMA(bars, i, longEmaWindow);
        double stdDev = computeStdDev(bars, i, longEmaWindow);
        double atr = computeEMA_ATR(bars, i, atrWindow);

        int signal = momentumSignal(currentPrice, shortEMA, mediumEMA);
        double z = zScore(currentPrice, longEMA, stdDev);
        double reactMult = reactionMultiplier(currentPrice, longEMA, stdDev);
        double size = reactMult * (riskUnit / (std::max(atr, 0.5)));

        if (!inTrade) {
            if (signal == 1 && z < 1.5) {
                inTrade = true;
                tradeDirection = 1;
                entryPrice = currentPrice;
                positionSize = size;
                entryIndex = i;
            } else if (signal == -1 && z > -1.5) {
                inTrade = true;
                tradeDirection = -1;
                entryPrice = currentPrice;
                positionSize = size;
                entryIndex = i;
            }
        } else {
            bool exitTrade = false;
            if (tradeDirection == 1) {
                if (z > 2.0 || signal != 1) {
                    exitTrade = true;
                }
            } else if (tradeDirection == -1) {
                if (z < -2.0 || signal != -1) {
                    exitTrade = true;
                }
            }

            if (exitTrade) {
                double tradePnL = 0.0;
                if (tradeDirection == 1) {
                    tradePnL = (currentPrice - entryPrice) * positionSize;
                } else if (tradeDirection == -1) {
                    tradePnL = (entryPrice - currentPrice) * positionSize;
                }

                // Record trade
                Trade trade;
                trade.entryTime = bars[entryIndex].timestamp;
                trade.exitTime = bars[i].timestamp;
                trade.entryPrice = entryPrice;
                trade.exitPrice = currentPrice;
                trade.direction = tradeDirection;
                trade.size = positionSize;
                trade.profit = tradePnL;
                trade.holdingPeriod = i - entryIndex;
                result.trades.push_back(trade);

                portfolio += tradePnL;
                inTrade = false;
                tradeDirection = 0;
                entryPrice = 0.0;
                positionSize = 0.0;
            }
        }

        double currentPortfolio = portfolio;
        if (inTrade) {
            if (tradeDirection == 1) {
                currentPortfolio += (currentPrice - entryPrice) * positionSize;
            } else if (tradeDirection == -1) {
                currentPortfolio += (entryPrice - currentPrice) * positionSize;
            }
        }

        result.timestamps.push_back(bars[i].timestamp);
        result.portfolioValues.push_back(currentPortfolio);
        if (result.portfolioValues.size() > 1) {
            double prev = result.portfolioValues[result.portfolioValues.size()-2];
            double ret = (currentPortfolio - prev) / prev;
            result.returns.push_back(ret);
        } else {
            result.returns.push_back(0.0);
        }
    }

    return result;
}

PerformanceMetrics computeMetrics(const BacktestResult &btResult, double initialCapital, int periodsPerYear) {
    PerformanceMetrics metrics = {};  // Initialize all members to 0
    const auto &values = btResult.portfolioValues;
    const auto &rets = btResult.returns;

    if (values.empty() || rets.empty()) {
        std::cerr << "No backtest data available for metrics." << std::endl;
        return metrics;
    }

    // Calculate return metrics
    double finalCapital = values.back();
    int totalPeriods = rets.size();
    metrics.annualizedReturn = std::pow(finalCapital / initialCapital, double(periodsPerYear) / totalPeriods) - 1;

    double sum = 0.0;
    for (double r : rets) { sum += r; }
    double mean = sum / rets.size();

    double sq_sum = 0.0;
    double downside_sq_sum = 0.0;
    for (double r : rets) {
        double deviation = r - mean;
        sq_sum += deviation * deviation;
        if (r < 0) {  // Only consider negative returns for Sortino
            downside_sq_sum += r * r;
        }
    }

    double stdev = (rets.size() > 1) ? std::sqrt(sq_sum / (rets.size()-1)) : 0.0;
    double downside_stdev = (rets.size() > 1) ? std::sqrt(downside_sq_sum / (rets.size()-1)) : 0.0;

    metrics.annualizedVolatility = stdev * std::sqrt(periodsPerYear);
    metrics.sharpeRatio = (stdev != 0) ? (mean / stdev) * std::sqrt(periodsPerYear) : 0.0;
    metrics.sortinoRatio = (downside_stdev != 0) ? (mean / downside_stdev) * std::sqrt(periodsPerYear) : 0.0;

    // Calculate drawdown
    double peak = values.front();
    double maxDD = 0.0;
    for (double v : values) {
        if (v > peak) peak = v;
        double drawdown = (peak - v) / peak;
        if (drawdown > maxDD) maxDD = drawdown;
    }
    metrics.maxDrawdown = maxDD;
    metrics.calmarRatio = (maxDD != 0) ? metrics.annualizedReturn / maxDD : 0.0;

    // Calculate trade metrics
    metrics.totalTrades = btResult.trades.size();

    double grossProfit = 0.0;
    double grossLoss = 0.0;
    int winCount = 0;
    int currentConsecutiveLosses = 0;
    double totalHoldingPeriods = 0.0;

    for (const auto& trade : btResult.trades) {
        totalHoldingPeriods += trade.holdingPeriod;

        if (trade.profit > 0) {
            winCount++;
            grossProfit += trade.profit;
            metrics.largestWin = std::max(metrics.largestWin, trade.profit);
            currentConsecutiveLosses = 0;
        } else {
            grossLoss -= trade.profit;  // Make loss positive for calculation
            metrics.largestLoss = std::min(metrics.largestLoss, trade.profit);
            currentConsecutiveLosses++;
            metrics.maxConsecutiveLosses = std::max(metrics.maxConsecutiveLosses,
                                                  (double)currentConsecutiveLosses);
        }
    }

    metrics.winningTrades = winCount;
    metrics.losingTrades = metrics.totalTrades - winCount;
    metrics.winRate = metrics.totalTrades > 0 ? (double)winCount / metrics.totalTrades : 0.0;
    metrics.averageWin = winCount > 0 ? grossProfit / winCount : 0.0;
    metrics.averageLoss = metrics.losingTrades > 0 ? grossLoss / metrics.losingTrades : 0.0;
    metrics.profitFactor = grossLoss > 0 ? grossProfit / grossLoss : 0.0;
    metrics.averageHoldingPeriod = metrics.totalTrades > 0 ? totalHoldingPeriods / metrics.totalTrades : 0.0;

    return metrics;
}

// -------------------------------
// CSV Output: Write Backtest Results to CSV
// -------------------------------
void outputCSV(const BacktestResult &btResult, const std::string &filename, const std::vector<Bar> &bars, int startIndex) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening CSV file for writing: " << filename << std::endl;
        return;
    }

    file << "Timestamp,Close,PortfolioValue,TradeDirection,PositionSize\n";
    for (size_t i = 0; i < btResult.timestamps.size(); i++) {
        file << btResult.timestamps[i] << ",";
        int barIndex = startIndex + i;
        if (barIndex < (int)bars.size()) {
            file << bars[barIndex].close << ",";
        } else {
            file << ",";
        }
        file << btResult.portfolioValues[i];

        // Find if there was a trade at this timestamp
        auto trade_it = std::find_if(btResult.trades.begin(), btResult.trades.end(),
            [&](const Trade& trade) {
                return trade.entryTime == btResult.timestamps[i] ||
                       trade.exitTime == btResult.timestamps[i];
            });

        if (trade_it != btResult.trades.end()) {
            file << "," << trade_it->direction << "," << trade_it->size;
        } else {
            file << ",,";
        }
        file << "\n";
    }

    // Write trade summary to a separate file
    std::string tradeSummaryFile = filename.substr(0, filename.find_last_of('.')) + "_trades.csv";
    std::ofstream tradeFile(tradeSummaryFile);
    if (tradeFile.is_open()) {
        tradeFile << "EntryTime,ExitTime,Direction,EntryPrice,ExitPrice,Size,Profit,HoldingPeriod\n";
        for (const auto& trade : btResult.trades) {
            tradeFile << trade.entryTime << ","
                     << trade.exitTime << ","
                     << trade.direction << ","
                     << trade.entryPrice << ","
                     << trade.exitPrice << ","
                     << trade.size << ","
                     << trade.profit << ","
                     << trade.holdingPeriod << "\n";
        }
        tradeFile.close();
        std::cout << "Trade summary written to " << tradeSummaryFile << std::endl;
    }

    file.close();
    std::cout << "CSV output written to " << filename << std::endl;
}

// -------------------------------
// Main: Process All Major Commodity Futures
// -------------------------------
int main() {
    std::vector<std::string> commoditySymbols = {
        // Agricultural
        "KEUSX",  // Wheat
        "ZSUSX",  // Soybeans
        "SBUSX", // Sugar
        "CTUSX", //Cotton

        // Metals
        "GCUSD",  // Gold
        "ALIUSD", // Aluminium

        // Energy
        "CLUSD",  // Crude Oil
        "NGUSD"   // Natural Gas

    };

    std::string timeframe = "30min";
    std::string fromDate = "2014-02-10";
    std::string toDate = "2024-03-10";
    std::string apikey = "hI1y91X0BXXmDIq7MVQWUGOXwyYufJp6";

    double initialCapital = 100000.0;
    int periodsPerYear = 3276;  // For 30-minute data: (6.5 hours * 2 periods/hour * 252 trading days)

    // Process each commodity symbol
    for (const auto &symbol : commoditySymbols) {
        std::cout << "\nProcessing " << symbol << "..." << std::endl;

        std::vector<Bar> bars = fetchHistoricalData(symbol, timeframe, fromDate, toDate, apikey);
        if (bars.empty()) {
            std::cerr << "No data retrieved for " << symbol << ". Skipping." << std::endl;
            continue;
        }
        std::cout << "Retrieved " << bars.size() << " bars of data." << std::endl;

        // Run backtest
        const int requiredBars = 100;
        if (bars.size() < (size_t)requiredBars) {
            std::cerr << "Not enough data for backtest. Skipping." << std::endl;
            continue;
        }

        BacktestResult btResult = runBacktest(bars, initialCapital);
        PerformanceMetrics metrics = computeMetrics(btResult, initialCapital, periodsPerYear);

        // Output results
        std::string filename = symbol + "_" + timeframe + "_backtest_results.csv";
        outputCSV(btResult, filename, bars, requiredBars);

        // Print performance metrics
        std::cout << "\n=== Performance Metrics for " << symbol << " ===" << std::endl;
        std::cout << "Return Metrics:" << std::endl;
        std::cout << "Annualized Return: " << (metrics.annualizedReturn * 100) << "%" << std::endl;
        std::cout << "Annualized Volatility: " << (metrics.annualizedVolatility * 100) << "%" << std::endl;
        std::cout << "Sharpe Ratio: " << metrics.sharpeRatio << std::endl;
        std::cout << "Sortino Ratio: " << metrics.sortinoRatio << std::endl;
        std::cout << "Max Drawdown: " << (metrics.maxDrawdown * 100) << "%" << std::endl;
        std::cout << "Calmar Ratio: " << metrics.calmarRatio << std::endl;

        std::cout << "\nTrade Metrics:" << std::endl;
        std::cout << "Total Trades: " << metrics.totalTrades << std::endl;
        std::cout << "Win Rate: " << (metrics.winRate * 100) << "%" << std::endl;
        std::cout << "Profit Factor: " << metrics.profitFactor << std::endl;
        std::cout << "Average Win: $" << metrics.averageWin << std::endl;
        std::cout << "Average Loss: $" << metrics.averageLoss << std::endl;
        std::cout << "Largest Win: $" << metrics.largestWin << std::endl;
        std::cout << "Largest Loss: $" << metrics.largestLoss << std::endl;
        std::cout << "Max Consecutive Losses: " << metrics.maxConsecutiveLosses << std::endl;
        std::cout << "Average Holding Period: " << metrics.averageHoldingPeriod << " bars" << std::endl;
    }

    return 0;
}
