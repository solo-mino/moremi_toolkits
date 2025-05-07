# In-Code Documentation Guide for `moremi_toolkits`

## 1. Purpose

Clear, concise, and consistent in-code documentation (docstrings and comments) is crucial for the maintainability, readability, and collaborative success of the `moremi_toolkits` project. This guide outlines the standards and best practices for documenting our Python code.

All public modules, classes, functions, and methods **must** have docstrings.

## 2. Docstring Style

We will adhere to the **Google Python Style Guide** for docstrings. This style is readable, comprehensive, and well-supported by documentation generation tools like Sphinx (with Napoleon extension).

**Key Features of Google Style Docstrings:**

*   Starts with a concise one-line summary.
*   Followed by a more detailed explanation if needed (after a blank line).
*   Uses specific sections like `Args:`, `Returns:`, `Raises:`, `Yields:`, `Attributes:`, `Examples:`.

### 2.1. General Docstring Format

```python
def my_function(arg1, arg2, arg3=None):
    """One-line summary of the function's purpose.

    More detailed explanation of what the function does, its behavior,
    and any important considerations for its use. Can span multiple lines.

    Args:
        arg1 (str): Description of the first argument.
        arg2 (int): Description of the second argument.
            Multi-line descriptions are indented.
        arg3 (Optional[bool]): Description of the optional third argument.
            Defaults to None.

    Returns:
        ReturnType: Description of the value returned. For example,
            a list of processed items.
        None: If the function does not explicitly return a value.

    Raises:
        ValueError: If `arg1` is an empty string.
        TypeError: If `arg2` is not an integer.

    Examples:
        >>> my_function("hello", 10)
        "Processed: hello with 10"
        >>> my_function("world", 5, arg3=True)
        "Processed: world with 5 (flagged)"

    Note:
        Any additional notes or important information.
    """
    # Function implementation
    if not arg1:
        raise ValueError("arg1 cannot be empty")
    if not isinstance(arg2, int):
        raise TypeError("arg2 must be an integer")
    # ...
    return f"Processed: {arg1} with {arg2}"
```

### 2.2. Module Docstrings

Every .py file should begin with a module-level docstring.

```python
"""One-line summary of the module's purpose.

This module provides utility functions for [specific domain, e.g., SMILES validation].
It includes functions for X, Y, and Z. It is used by [other modules/components].
"""

# Imports follow the module docstring
import os
# ...
```
### 2.3. Class Docstrings
```python
Classes should have a docstring immediately below the class definition.
Include an Attributes: section for public attributes if they are not self-evident from __init__ arguments.

class MyClass:
    """One-line summary of the class's purpose.

    More detailed explanation of the class, its responsibilities,
    and how it should be instantiated and used.

    Attributes:
        public_attr1 (str): Description of a public attribute.
        _internal_attr (int): Internal attributes (prefixed with underscore)
            should generally not be in the public docstring unless necessary
            for understanding subclassing.
    """

    def __init__(self, param1, param2):
        """Initializes MyClass.

        Args:
            param1 (str): Description of parameter 1.
            param2 (list[int]): Description of parameter 2.
        """
        self.public_attr1 = param1
        self._internal_attr = len(param2)
        # ...

    def public_method(self, arg):
        """Summary of what the public method does.

        Args:
            arg (Any): Description of the argument.

        Returns:
            bool: True if successful, False otherwise.
        """
        # ...
        return True

    def _private_method(self):
        """Docstrings for private methods (single underscore) are encouraged
        for maintainer clarity but are not part of the public API.
        """
        # ...
```

### 2.4. Type Hinting

- Always use type hints for function arguments, return values, and class attributes where feasible.

- Type hints complement docstrings by providing machine-readable type information. The docstring then explains the meaning and constraints of those types.

- Use types from the typing module (e.g., Optional, List, Dict, Tuple, Callable, Any).

```python
from typing import List, Optional, Dict

def process_data(data: List[Dict[str, str]], threshold: Optional[float] = None) -> int:
    """Processes a list of data dictionaries.

    Args:
        data: A list of dictionaries, where each dictionary represents an item.
        threshold: An optional float value to filter items. Defaults to None.

    Returns:
        The count of items processed.
    """
    # ...
    return 0
```

## 3. Inline Comments
Inline comments (#) should be used to explain complex, non-obvious, or tricky parts of the code.

- **Clarity over Quantity**: Don't comment on obvious code. Comments should add understanding that the code itself cannot easily convey.

```python
# BAD: Obvious comment
x = x + 1 # Increment x

# GOOD: Explains the "why" or a complex condition
# We need to adjust the index here because the external API is 1-based.
adjusted_index = original_index - 1
```
- Keep Comments Up-to-Date: If you change the code, ensure any related comments are also updated or removed. Stale comments are worse than no comments.

- `TODO` and `FIXME` Comments:

    - Use `TODO`: to mark areas that need future work or improvement. Include your name/initials and date if helpful, or a ticket/issue number.
        - `TODO` (J.Doe 2024-01-15): Refactor this to use the new foobar_service.
        - `TODO`: Add more robust error handling here (see issue #123).

    - Use `FIXME`: to mark known bugs or areas that need correction.
        - `FIXME`: This calculation is incorrect under condition X.

- Placement:
    - Comments explaining a block of code should typically precede it.
    - Short comments on a single line can follow the statement.

## 4. Writing Style for Documentation

- **Clear and Concise**: Use simple, direct language. Avoid jargon where possible, or explain it if necessary.

- **Audience**: Write for other developers (including your future self) who may not have the same context you do.

- **Imperative Mood**: For function/method summaries, use the imperative mood (e.g., "Return the sum..." not "Returns the sum..." or "This function returns the sum...").

- **Complete Sentences**: Use complete sentences and proper grammar/punctuation.

- **Consistency**: Follow the chosen style (Google) consistently throughout the project.

## 5. Reviewing Documentation

Documentation should be treated as an integral part of the code. During code reviews:

- Check for the presence and correctness of docstrings for all public APIs.

- Verify that docstrings accurately reflect the code's behavior.

- Ensure type hints are present and correct.

- Assess the clarity and usefulness of inline comments.