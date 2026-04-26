"""
Desktop GUI wrapper for fullPlasmidSeq_demix_RFW.py.

Run directly with `python gui.py`, or build a standalone bundle with build_app.sh.
"""
import os
import sys
import queue
import logging
import threading
import traceback
from pathlib import Path
from tkinter import filedialog, messagebox, PhotoImage

import customtkinter as ctk

from setup_tools import ensure_tools_available, ToolsMissingError
import fullPlasmidSeq_demix_RFW as demix


def _resource_path(relative: str) -> Path:
    """Resolve a bundled resource both in dev and inside a PyInstaller bundle."""
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    return base / relative


ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")


class QueueLogHandler(logging.Handler):
    """Push log records to a queue for the UI thread to drain."""
    def __init__(self, log_queue: queue.Queue):
        super().__init__()
        self.log_queue = log_queue

    def emit(self, record):
        try:
            self.log_queue.put(self.format(record))
        except Exception:
            self.handleError(record)


class DemixApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("Plasmid Demix")
        self.geometry("780x640")
        self.minsize(680, 560)
        self._set_window_icon()

        self.log_queue: queue.Queue[str] = queue.Queue()
        self.worker: threading.Thread | None = None
        self._cancelling = False

        self._build_ui()
        self.after(100, self._drain_log_queue)

    def _set_window_icon(self):
        """Set the title-bar/dock icon if assets/icon.png exists."""
        icon_png = _resource_path("assets/icon.png")
        if icon_png.exists():
            try:
                self._icon_image = PhotoImage(file=str(icon_png))
                self.iconphoto(True, self._icon_image)
            except Exception:
                pass

    def _build_ui(self):
        pad = {"padx": 12, "pady": 6}

        header = ctk.CTkLabel(
            self, text="Plasmid Demix",
            font=ctk.CTkFont(size=22, weight="bold"),
        )
        header.pack(pady=(16, 4))

        subtitle = ctk.CTkLabel(
            self,
            text="Demix pooled nanopore plasmid reads and build per-plasmid consensus sequences.",
            font=ctk.CTkFont(size=12),
            text_color=("gray30", "gray70"),
        )
        subtitle.pack(pady=(0, 12))

        form = ctk.CTkFrame(self)
        form.pack(fill="x", **pad)
        form.grid_columnconfigure(1, weight=1)

        self.excel_var = ctk.StringVar()
        self.fastq_var = ctk.StringVar()
        self.ref_var = ctk.StringVar()
        self.out_var = ctk.StringVar()

        self._add_path_row(form, 0, "Excel file:", self.excel_var, self._pick_excel)
        self._add_path_row(form, 1, "FASTQ folder:", self.fastq_var, lambda: self._pick_dir(self.fastq_var))
        self._add_path_row(form, 2, "Reference folder (optional):", self.ref_var, lambda: self._pick_dir(self.ref_var))
        self._add_path_row(form, 3, "Output folder:", self.out_var, lambda: self._pick_dir(self.out_var))

        opts = ctk.CTkFrame(self)
        opts.pack(fill="x", **pad)
        opts.grid_columnconfigure(3, weight=1)

        self.quick_var = ctk.BooleanVar(value=False)
        self.keep_temp_var = ctk.BooleanVar(value=False)
        self.threads_var = ctk.IntVar(value=4)

        ctk.CTkCheckBox(opts, text="Quick method (skip de novo assembly)",
                        variable=self.quick_var).grid(row=0, column=0, columnspan=2, sticky="w", padx=10, pady=8)
        ctk.CTkCheckBox(opts, text="Keep temp files",
                        variable=self.keep_temp_var).grid(row=0, column=2, sticky="w", padx=10, pady=8)

        ctk.CTkLabel(opts, text="Threads:").grid(row=1, column=0, sticky="w", padx=10, pady=8)
        threads_spin = ctk.CTkOptionMenu(
            opts,
            values=[str(n) for n in (1, 2, 4, 6, 8, 12, 16)],
            variable=ctk.StringVar(value="4"),
            command=lambda v: self.threads_var.set(int(v)),
            width=80,
        )
        threads_spin.grid(row=1, column=1, sticky="w", padx=4, pady=8)

        action = ctk.CTkFrame(self, fg_color="transparent")
        action.pack(fill="x", **pad)
        self.run_button = ctk.CTkButton(action, text="Run", width=120, command=self._on_run)
        self.run_button.pack(side="right", padx=4)
        self.stop_button = ctk.CTkButton(
            action, text="Stop", width=120, command=self._on_stop,
            fg_color="#b3261e", hover_color="#8c1d18", state="disabled",
        )
        self.stop_button.pack(side="right", padx=4)
        self.status_label = ctk.CTkLabel(action, text="Idle", anchor="w")
        self.status_label.pack(side="left", padx=4)

        log_frame = ctk.CTkFrame(self)
        log_frame.pack(fill="both", expand=True, **pad)
        ctk.CTkLabel(log_frame, text="Log", anchor="w").pack(fill="x", padx=8, pady=(6, 0))
        self.log_box = ctk.CTkTextbox(log_frame, wrap="none", font=ctk.CTkFont(family="Menlo", size=11))
        self.log_box.pack(fill="both", expand=True, padx=8, pady=8)
        self.log_box.configure(state="disabled")

    def _add_path_row(self, parent, row, label, var, picker):
        ctk.CTkLabel(parent, text=label, width=180, anchor="w").grid(row=row, column=0, sticky="w", padx=10, pady=6)
        entry = ctk.CTkEntry(parent, textvariable=var)
        entry.grid(row=row, column=1, sticky="ew", padx=4, pady=6)
        ctk.CTkButton(parent, text="Browse…", width=90, command=picker).grid(row=row, column=2, padx=10, pady=6)

    def _pick_excel(self):
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if path:
            self.excel_var.set(path)

    def _pick_dir(self, var: ctk.StringVar):
        path = filedialog.askdirectory(title="Select folder")
        if path:
            var.set(path)

    def _validate(self) -> tuple[bool, str]:
        if not self.excel_var.get() or not Path(self.excel_var.get()).is_file():
            return False, "Excel file is required."
        if not self.fastq_var.get() or not Path(self.fastq_var.get()).is_dir():
            return False, "FASTQ folder is required."
        if not self.out_var.get():
            return False, "Output folder is required."
        if self.ref_var.get() and not Path(self.ref_var.get()).is_dir():
            return False, "Reference folder path is invalid."
        return True, ""

    def _on_run(self):
        if self.worker and self.worker.is_alive():
            return
        ok, msg = self._validate()
        if not ok:
            messagebox.showerror("Missing input", msg)
            return

        self._clear_log()
        self._cancelling = False
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal", text="Stop")
        self.status_label.configure(text="Checking bioinformatics tools…")
        self.worker = threading.Thread(target=self._run_pipeline, daemon=True)
        self.worker.start()

    def _on_stop(self):
        if not (self.worker and self.worker.is_alive()):
            return
        if self._cancelling:
            return
        self._cancelling = True
        self.stop_button.configure(state="disabled", text="Stopping…")
        self._set_status("Cancelling — terminating running tools…")
        self.log_queue.put("[stop] Cancellation requested by user.")
        threading.Thread(target=demix.request_cancel, daemon=True).start()

    def _run_pipeline(self):
        handler = QueueLogHandler(self.log_queue)
        handler.setFormatter(logging.Formatter("%(asctime)s — %(levelname)s — %(message)s"))
        root_logger = logging.getLogger()
        root_logger.addHandler(handler)
        try:
            try:
                ensure_tools_available(log=lambda m: self.log_queue.put(m))
            except ToolsMissingError as e:
                self.log_queue.put(f"[setup] {e}")
                self._set_status(f"Tools missing: {e}")
                return

            self._set_status("Running…")
            output_dir = self.out_var.get()
            os.makedirs(output_dir, exist_ok=True)

            demix.main(
                excel_file=self.excel_var.get(),
                ref_dir=self.ref_var.get() or "",
                fastq_dir=self.fastq_var.get(),
                output_dir=output_dir,
                keep_temp=self.keep_temp_var.get(),
                threads=int(self.threads_var.get()),
                quick_method=self.quick_var.get(),
            )
            self._set_status("Done.")
            self.log_queue.put("[done] Pipeline finished successfully.")
        except demix.PipelineCancelled:
            self._set_status("Cancelled.")
            self.log_queue.put("[stop] Pipeline cancelled.")
        except Exception as e:
            self.log_queue.put("[error] " + "".join(traceback.format_exception_only(type(e), e)).strip())
            self.log_queue.put(traceback.format_exc())
            self._set_status(f"Failed: {e}")
        finally:
            root_logger.removeHandler(handler)
            self.after(0, self._reset_buttons)

    def _reset_buttons(self):
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled", text="Stop")
        self._cancelling = False

    def _set_status(self, text: str):
        self.after(0, lambda: self.status_label.configure(text=text))

    def _clear_log(self):
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")

    def _drain_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self.log_box.configure(state="normal")
                self.log_box.insert("end", msg + "\n")
                self.log_box.see("end")
                self.log_box.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(120, self._drain_log_queue)


def main():
    app = DemixApp()
    app.mainloop()


if __name__ == "__main__":
    main()
