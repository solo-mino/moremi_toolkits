"""
Internal Utilities for Connector Modules.

This module contains helper functions and classes that are shared among
the different connector modules within the `moremi_biokit.connectors`
sub-package. These utilities are not intended for direct public use
outside of this sub-package.

Functionalities might include:
- Common HTTP request wrappers (e.g., with consistent error handling, retries, timeout settings).
- Helper functions for parsing common API response formats (e.g., XML, specific JSON structures).
- Rate limiting helpers or decorators.
- User-agent string definition.

The leading underscore in the filename (`_utils.py`) indicates that this
is an internal implementation detail.
"""

import time
import requests
from typing import Optional, Dict, Any, Union

# Default Moremi Biokit User-Agent
MOREMI_BIOKIT_USER_AGENT = "Moremi-Biokit/0.1 (https://github.com/solo-mino/moremi_toolkits)"

class APIRequestError(Exception):
    """Custom exception for API request failures after retries."""
    def __init__(self, message: str, status_code: Optional[int] = None, url: Optional[str] = None):
        super().__init__(message)
        self.status_code = status_code
        self.url = url

    def __str__(self) -> str:
        return f"APIRequestError: {super().__str__()} (Status: {self.status_code}, URL: {self.url})"


def make_api_request(
    url: str,
    method: str = "GET",
    params: Optional[Dict[str, Any]] = None,
    data: Optional[Union[Dict[str, Any], str]] = None,
    json_data: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = 30, # seconds
    retries: int = 3,
    delay: int = 5, # seconds between retries
    stream: bool = False
) -> requests.Response:
    """Makes an HTTP request with error handling and retries.

    Args:
        url (str): The URL for the request.
        method (str, optional): HTTP method (GET, POST, etc.). Defaults to "GET".
        params (Optional[Dict[str, Any]], optional): URL parameters. Defaults to None.
        data (Optional[Union[Dict[str, Any], str]], optional): Data to send in the body (for POST/PUT). Defaults to None.
        json_data (Optional[Dict[str, Any]], optional): JSON data to send in the body. Defaults to None.
        headers (Optional[Dict[str, str]], optional): HTTP headers. Defaults to None.
                                                   A default User-Agent is added if not provided.
        timeout (int, optional): Request timeout in seconds. Defaults to 30.
        retries (int, optional): Number of retries for transient errors. Defaults to 3.
        delay (int, optional): Delay in seconds between retries. Defaults to 5.
        stream (bool, optional): If True, the response content will not be immediately downloaded.
                                 Defaults to False.

    Returns:
        requests.Response: The response object.

    Raises:
        APIRequestError: If the request fails after all retries, or for non-transient HTTP errors.
    """
    effective_headers = {"User-Agent": MOREMI_BIOKIT_USER_AGENT}
    if headers:
        effective_headers.update(headers)

    for attempt in range(retries + 1):
        try:
            response = requests.request(
                method=method.upper(),
                url=url,
                params=params,
                data=data,
                json=json_data,
                headers=effective_headers,
                timeout=timeout,
                stream=stream
            )
            # Raise HTTPError for bad responses (4xx or 5xx)
            # We handle 5xx as potentially transient below, but 4xx (client errors) are usually not.
            # For 4xx errors, we might not want to retry.
            if 400 <= response.status_code < 500 and response.status_code not in [408, 429]: # 408: Timeout, 429: Too Many Requests
                response.raise_for_status() # Will raise HTTPError for these specific client errors

            # For server errors (5xx) or specific retryable client errors (408, 429), we attempt retries.
            if response.status_code >= 500 or response.status_code in [408, 429]:
                response.raise_for_status() # This will be caught by requests.exceptions.HTTPError for retry logic

            return response  # Success
        
        except requests.exceptions.HTTPError as e:
            # If it's a client error that we decided not to retry (e.g. 404 Not Found), raise immediately.
            if 400 <= e.response.status_code < 500 and e.response.status_code not in [408, 429]:
                raise APIRequestError(
                    f"Client error: {e.response.status_code} {e.response.reason} for URL: {url}",
                    status_code=e.response.status_code,
                    url=url
                ) from e
            
            # For server errors or retryable client errors
            print(f"Warning: HTTP error {e.response.status_code} for {url}. Retrying ({attempt + 1}/{retries + 1})...")
            if attempt >= retries:
                raise APIRequestError(
                    f"HTTP error after {retries} retries: {e.response.status_code} {e.response.reason}",
                    status_code=e.response.status_code,
                    url=url
                ) from e
        
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            print(f"Warning: Request failed ({type(e).__name__}) for {url}. Retrying ({attempt + 1}/{retries + 1})...")
            if attempt >= retries:
                raise APIRequestError(
                    f"Request failed after {retries} retries: {type(e).__name__}",
                    url=url
                ) from e
        
        except requests.exceptions.RequestException as e: # Catch any other request-related errors
            print(f"Warning: Unexpected request error ({type(e).__name__}) for {url}. Retrying ({attempt + 1}/{retries + 1})...")
            if attempt >= retries:
                raise APIRequestError(
                    f"Unexpected request error after {retries} retries: {type(e).__name__}",
                    url=url
                ) from e
        
        if attempt < retries:
            time.sleep(delay * (attempt + 1)) # Exponential backoff can be added here

    # Should not be reached if logic is correct, but as a fallback:
    raise APIRequestError("Request failed after all retries due to an unknown issue.", url=url)

# Function definitions like:
# def _make_api_request(url: str, params: Optional[Dict] = None, method: str = "GET", ...) -> requests.Response:
# def _handle_api_error(response: requests.Response):
# ...