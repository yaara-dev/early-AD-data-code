@echo off
REM Simple script to start a local web server for Windows
REM This allows the HTML to automatically load CSV files from the data/csvs folder

echo ==========================================
echo   ExSeq Brain AD Visualization Server
echo ==========================================
echo.
echo Starting local web server...
echo.

REM Open browser after 2 seconds
start /b timeout /t 2 /nobreak >nul && start http://localhost:8000

echo The visualization will be available at:
echo   http://localhost:8000
echo.
echo Your browser should open automatically in a few seconds.
echo Press Ctrl+C to stop the server
echo.
echo ==========================================
echo.

REM Try Python 3 first
python -m http.server 8000 2>nul
if %errorlevel% neq 0 (
    REM Try Python 2
    python -m SimpleHTTPServer 8000 2>nul
    if %errorlevel% neq 0 (
        echo Error: Could not find Python to start a web server.
        echo Please install Python 3 from https://www.python.org/
        echo Or use the file input in the HTML to load CSV files manually.
        echo.
        pause
        exit /b 1
    )
)

pause

