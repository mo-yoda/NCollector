"""
Tests for dialogs.py — ask_user_parameter modal dialog.

Uses real tkinter with root.after() to simulate user input on the event loop.
These tests will be skipped in headless CI environments (no display).
"""
import pytest

try:
    import tkinter as tk
    # Test if display is available
    _test_root = tk.Tk()
    _test_root.destroy()
    TK_AVAILABLE = True
except (tk.TclError, Exception):
    TK_AVAILABLE = False

pytestmark = pytest.mark.skipif(not TK_AVAILABLE, reason="No display available for tkinter")

from dialogs import ask_user_parameter


@pytest.fixture
def root():
    """Create and destroy a tkinter root window for each test."""
    r = tk.Tk()
    r.withdraw()  # Hide the root window
    yield r
    r.destroy()


def _simulate_input_and_confirm(root, dialog_title, value):
    """
    Schedules typing a value into the dialog's entry and clicking Confirm.
    Must be called BEFORE ask_user_parameter (which blocks on wait_window).
    """
    def _action():
        # Find the Toplevel dialog by title
        for widget in root.winfo_children():
            if isinstance(widget, tk.Toplevel) and widget.title() == dialog_title:
                # Find the Entry widget
                for child in widget.winfo_children():
                    if isinstance(child, tk.Entry):
                        child.delete(0, tk.END)
                        child.insert(0, str(value))
                        break
                # Find and click Confirm button (first button in the button frame)
                for child in widget.winfo_children():
                    if isinstance(child, tk.Frame):
                        for btn in child.winfo_children():
                            if isinstance(btn, tk.Button) and btn.cget("text") == "Confirm":
                                btn.invoke()
                                return
    root.after(100, _action)


def _simulate_skip(root, dialog_title):
    """Schedules clicking the Skip button."""
    def _action():
        for widget in root.winfo_children():
            if isinstance(widget, tk.Toplevel) and widget.title() == dialog_title:
                for child in widget.winfo_children():
                    if isinstance(child, tk.Frame):
                        for btn in child.winfo_children():
                            if isinstance(btn, tk.Button) and btn.cget("text") == "Skip":
                                btn.invoke()
                                return
    root.after(100, _action)


def _simulate_close(root, dialog_title):
    """Schedules closing the dialog via window manager (X button)."""
    def _action():
        for widget in root.winfo_children():
            if isinstance(widget, tk.Toplevel) and widget.title() == dialog_title:
                widget.destroy()
                return
    root.after(100, _action)


# ---------------------------------------------------------------------------
# Tests: Integer input
# ---------------------------------------------------------------------------

class TestIntInput:

    def test_valid_int(self, root):
        """Valid integer input should return int."""
        title = "Test Int"
        _simulate_input_and_confirm(root, title, "42")
        result = ask_user_parameter(root, title=title, message="Enter int", input_type="int")
        assert result == 42
        assert isinstance(result, int)

    def test_negative_int(self, root):
        """Negative integer should work."""
        title = "Test Neg Int"
        _simulate_input_and_confirm(root, title, "-5")
        result = ask_user_parameter(root, title=title, message="Enter int", input_type="int")
        assert result == -5

    def test_zero(self, root):
        """Zero should be a valid int (not treated as empty)."""
        title = "Test Zero"
        _simulate_input_and_confirm(root, title, "0")
        result = ask_user_parameter(root, title=title, message="Enter int", input_type="int")
        assert result == 0


# ---------------------------------------------------------------------------
# Tests: Float input
# ---------------------------------------------------------------------------

class TestFloatInput:

    def test_valid_float(self, root):
        """Valid float input should return float."""
        title = "Test Float"
        _simulate_input_and_confirm(root, title, "3.14")
        result = ask_user_parameter(root, title=title, message="Enter float", input_type="float")
        assert result == pytest.approx(3.14)
        assert isinstance(result, float)

    def test_int_as_float(self, root):
        """Integer string should work for float input."""
        title = "Test Int Float"
        _simulate_input_and_confirm(root, title, "5")
        result = ask_user_parameter(root, title=title, message="Enter float", input_type="float")
        assert result == pytest.approx(5.0)

    def test_scientific_notation(self, root):
        """Scientific notation should parse as float."""
        title = "Test Sci"
        _simulate_input_and_confirm(root, title, "1e-3")
        result = ask_user_parameter(root, title=title, message="Enter float", input_type="float")
        assert result == pytest.approx(0.001)


# ---------------------------------------------------------------------------
# Tests: String input
# ---------------------------------------------------------------------------

class TestStrInput:

    def test_valid_string(self, root):
        """String input should return the raw string."""
        title = "Test Str"
        _simulate_input_and_confirm(root, title, "hello world")
        result = ask_user_parameter(root, title=title, message="Enter str", input_type="str")
        assert result == "hello world"

    def test_numeric_as_string(self, root):
        """Numeric input with type='str' should stay as string."""
        title = "Test Num Str"
        _simulate_input_and_confirm(root, title, "42")
        result = ask_user_parameter(root, title=title, message="Enter str", input_type="str")
        assert result == "42"
        assert isinstance(result, str)


# ---------------------------------------------------------------------------
# Tests: Skip / Cancel
# ---------------------------------------------------------------------------

class TestSkipAndCancel:

    def test_skip_returns_none(self, root):
        """Clicking Skip should return None."""
        title = "Test Skip"
        _simulate_skip(root, title)
        result = ask_user_parameter(root, title=title, message="Skip me", input_type="int")
        assert result is None

    def test_close_window_returns_none(self, root):
        """Closing the dialog window should return None."""
        title = "Test Close"
        _simulate_close(root, title)
        result = ask_user_parameter(root, title=title, message="Close me", input_type="int")
        assert result is None


# ---------------------------------------------------------------------------
# Tests: Default value
# ---------------------------------------------------------------------------

class TestDefaultValue:

    def test_default_prefilled(self, root):
        """Default value should be pre-filled; confirming it returns the default."""
        title = "Test Default"
        # Simulate confirming without changing the pre-filled value
        def _confirm_default():
            for widget in root.winfo_children():
                if isinstance(widget, tk.Toplevel) and widget.title() == title:
                    for child in widget.winfo_children():
                        if isinstance(child, tk.Frame):
                            for btn in child.winfo_children():
                                if isinstance(btn, tk.Button) and btn.cget("text") == "Confirm":
                                    btn.invoke()
                                    return
        root.after(100, _confirm_default)
        result = ask_user_parameter(root, title=title, message="Use default",
                                    input_type="int", default=7)
        assert result == 7

    def test_default_none(self, root):
        """No default should leave entry empty; entering a value should work."""
        title = "Test No Default"
        _simulate_input_and_confirm(root, title, "99")
        result = ask_user_parameter(root, title=title, message="No default",
                                    input_type="int", default=None)
        assert result == 99