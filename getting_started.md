# Getting Started

General instructions on the use of these tutorials

## Summary

This document describes how the directory structure should be set up, how commands can be copied and executed from tutorial instructions, how text files can be used and edited, and how files required for tutorials can be downloaded from the tutorial repository.

## Table of contents

* [Directory structure](#directory_structure)
* [Copying and pasting commands](#copy_pasting)
* [Working with text files](#text_files)
* [Downloading tutorial files](#downloading)

<a name="directory_structure"></a>
## Directory structure

Following the tutorials in this collection may be easiest if separate directories are used for each tutorial. These directories could be named after the tutorial, so for example, to follow tutorial [Substitution Model Selection](substitution_model_selection/README.md), a directory named `substitution_model_selection` could be generated. It should not matter where on your computer the directory is placed as long as you do have the rights to execute programs and commands in that directory. After the directory has been generated, let's say, in `/USERNAME/Desktop/`, it would be best to navigate to that directory both with a file manager program (such as Finder on Mac OS X, Konqueror on Linux, or Explorer on Windows) and with a console program (like the program called Terminal on Mac OS X and Linux or the Windows Console on Windows). The console program will be important to execute commands, which will sometimes be required during the tutorials. You will recognize such commands in the tutorials by monospace font, gray background, and an outline, like for example this command:

		pwd
		
To execute commands like the above, the "Enter" key always needs to be hit after writing the command. In general, steps that should be executed by you during the tutorial are marked with a bullet point, such as this:

* Unless you have already navigated into the tutorial directory with your console program (because you already know how to do that), do so now. If the directory would be named and placed as described above, the command to do so would be this one (assuming that BASH is installed and used; see [Requirements](requirements.md)):

		cd /USERNAME/Desktop/substitution_model_selection
		
	Because your actual user name is unlikely to be "USERNAME" and you may have placed the directory not on the Desktop but elsewhere, you will need to adjust the above command before executing it.

* To verify that you have successfully navigated to the tutorial directory, you could type this following command, which is an abbreviation for "print working directory":

		pwd
		
	The console should then print a line with the path and the name of the tutorial directory, such as "/USERNAME/Desktop/substitution\_model\_selection".

<a name="copy_pasting"></a>
## Copying and pasting commands

After you've navigated to the tutorial directory with the console program, it should be possible to execute all commands that are mentioned in the tutorial instructions exactly as they are specified. All one-line commands can be copied from the instructions and pasted into the console window, but care should be taken not to include whitespace symbols before the first letter or after the last letter of the command. For multi-line commands, however, it is safer to copy and paste each line individually, and hit the Enter key each time a line has been copied (again, whitespace before and after the first and last letters of the line should not be copied).
	
* Try to copy and paste this command into the console window, and then execute it:

		for n in {1..10}
		do
			echo ${n}
		done

	This should result in the numbers 1 to 10 being written line by line in the console window; if this is not the resulting output, something must have gone wrong with copy-pasting.

<a name="text_files"></a>
## Working with text files

In some cases, text files will need to be written or edited, and for this, a text editor of some sort will be required. For example, instead of copying many lines of commands into the console window one by one, these could also be copied into a text file, the file could be named, e.g., `script.sh`, and then all commands from that file could be executed jointly the command

		bash script.sh
		
To write and edit text files, either a command-line text editor or one with a graphical user interface could be used. There are many suitable options in the latter category. Programs like TextEdit (on Mac OS X) or Notepad (on Windows) could be used, but only if the default settings are changed so that plain text instead of "rich text" is written. Using Microsoft Word would guarantee trouble sooner rather than later. More convenient than any of these, however, are those that include syntax highlighting, such as [BBEdit](https://www.barebones.com/products/textwrangler/) for Mac OS X, [Notepad++](https://notepad-plus-plus.org) for Windows, or [Geany](https://www.geany.org) and [Atom](https://atom.io) that run on all three platforms.

If you want to avoid installing any of these programs, you could also use one of the text editors that are usually available within the console, like Vim or Emacs. No installation should be required for these tools, but their use may not be as intuitive as that of text editors with graphical user interfaces. Short tutorials are available online for [Vim](https://www.howtoforge.com/vim-basics) and [Emacs](https://www.digitalocean.com/community/tutorials/how-to-use-the-emacs-editor-in-linux).

* Use the text editor of your choice to write the following four lines:

		for n in {1..10}
		do
			echo ${n}
		done

* Save the text to a new file named `script.sh` that is placed in the tutorial directory.

* Execute the commands in file `script.sh` with this command:

		bash script.sh
		
	This should again result in the numbers 1 to 10 being written on separate lines.

<a name="downloading"></a>
## Downloading tutorial files

Links to input files and scripts are included in all tutorials and to run a given tutorial, the linked files should be downloaded and placed in the tutorial directory. Most often, these files are hosted in the same repository, and will not automatically download when the link is clicked (assuming that the tutorial instructions are opened in a web browser).

[`16s_filtered.nex`](substitution_model_selection/data/16s_filtered.nex)