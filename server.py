#!/usr/bin/env python3
"""
Simple HTTP server with no-cache headers for CSV files to prevent browser caching issues.
"""

import http.server
import socketserver
from pathlib import Path

class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Add no-cache headers for CSV files
        if self.path.endswith('.csv'):
            self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
            self.send_header('Pragma', 'no-cache')
            self.send_header('Expires', '0')
        super().end_headers()

if __name__ == "__main__":
    PORT = 8000
    
    with socketserver.TCPServer(("", PORT), NoCacheHTTPRequestHandler) as httpd:
        print(f"Server running at http://localhost:{PORT}/")
        print("CSV files will be served with no-cache headers")
        httpd.serve_forever()

