import tkinter as tk


def ask_user_parameter(parent,
                       title: str,
                       message: str,
                       input_type: str = "int",
                       default=None):
    """
    Reusable modal dialog for requesting a single parameter from the user.

    Args:
        parent: Tkinter parent window.
        title: Dialog window title.
        message: Prompt message displayed to user.
        input_type: Expected type — "int", "float", or "str".
        default: Default value pre-filled in entry field.

    Returns:
        The user's input (cast to input_type), or None if skipped/cancelled.
    """
    result = [None]  # Mutable container for closure

    dialog = tk.Toplevel(parent)
    dialog.title(title)
    dialog.geometry("450x220")
    dialog.transient(parent)
    dialog.grab_set()

    # Message
    tk.Label(dialog, text=message, wraplength=400, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15)

    # Entry
    entry_var = tk.StringVar(value=str(default) if default is not None else "")
    entry = tk.Entry(dialog, textvariable=entry_var, width=15, font=("Arial", 11))
    entry.pack(pady=5)
    entry.focus_set()

    # Validation feedback
    error_label = tk.Label(dialog, text="", fg="red", font=("Arial", 9))
    error_label.pack()

    def on_confirm():
        raw = entry_var.get().strip()
        if not raw:
            error_label.config(text="Please enter a value.")
            return
        try:
            if input_type == "int":
                result[0] = int(raw)
            elif input_type == "float":
                result[0] = float(raw)
            else:
                result[0] = raw
            dialog.destroy()
        except ValueError:
            error_label.config(text=f"Invalid input. Expected a {input_type} value.")

    def on_skip():
        dialog.destroy()

    # Buttons
    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=15)

    tk.Button(btn_frame, text="Confirm", command=on_confirm,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Skip", command=on_skip,
              padx=10).pack(side="left", padx=10)

    # Allow Enter key to confirm
    entry.bind("<Return>", lambda e: on_confirm())

    parent.wait_window(dialog)
    return result[0]