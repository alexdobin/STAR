# Contributing to STAR

The following is a set of guidelines for contributing to STAR hosted in the GitHub: https://github.com/alexdobin/STAR/. 
These are mostly guidelines, not rules. 
Use your best judgment, and feel free to propose changes to this document in a pull request.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[How Can I Contribute?](#how-can-i-contribute)
  * [Ask a question](#ask-a-question)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Pull Requests](#pull-requests)

## Code of Conduct

This project and everyone participating in it is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.

## How Can I Contribute?

### Ask a question

Please do not file an issue to ask a question.
The GitHub issue tracker is intended for bug reports and feature requests.
We have an official discussion forum where the community chimes in with helpful advice if you have questions.

* [STAR discussion forum](https://groups.google.com/forum/#!forum/rna-star)

### Reporting Bugs

This section guides you through submitting a bug report for STAR. 
Following these guidelines helps maintainers and the community understand your report, 
reproduce the behavior, and find related reports.
If you find a **Closed** issue that seems like it is the same thing that you're experiencing, 
open a new issue and include a link to the original issue in the body of your new one.

#### Before Submitting A Bug Report

* You might be able to find the cause of the problem and fix things yourself. 
Most importantly, check if you can reproduce the problem in the latest version of STAR
and if the problem happens when you run with mostly default parameters.
* Check the Log.out file for ERROR/WARNING/SOLUTION messages.
* Perform a through STAR GitHub issues to see if the problem has already been reported. 
If it has **and the issue is still open**, add a comment to the existing issue instead of opening a new one.

#### How Do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/) on https://github.com/alexdobin/STAR/issues.
Explain the problem and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible. For example, start by explaining how you run STAR, e.g. which command exactly you used in the terminal. When listing steps, **don't just say what you did, but explain how you did it.
* **Log.out**. Attach the Log.out file generated in the failed run. This file contains a lot of useful debugging information and is a starting point
* **Provide specific examples to demonstrate the steps**. Include links to files or GitHub projects, or copy/pasteable snippets, which you use in those examples. If you're providing snippets in the issue, use [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the behavior you observed after following the steps** and point out what exactly is the problem with that behavior.
* **Explain which behavior you expected to see instead and why.**
* **System information**. In many cases, the problems are associated with the hardware configurations. Provide a brief description of the CPU(s), RAM, storage. 
* **Can you reliably reproduce the issue?** If not, provide details about how often the problem happens and under which conditions it normally happens.

Include details about your configuration and environment:

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for Atom, including completely new features and minor improvements to existing functionality. Following these guidelines helps maintainers and the community understand your suggestion and find related suggestions.

#### Before Submitting An Enhancement Suggestion

* **Check the STAR manual and Release Notes** — you might discover that the enhancement is already available. Most importantly, check if you're using the latest version of STAR.

* **Perform a cursory search** to see if the enhancement has already been suggested. If it has, add a comment to the existing issue instead of opening a new one.

#### How Do I Submit A (Good) Enhancement Suggestion?

Feature requests and enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/) on https://github.com/alexdobin/STAR/issues.

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Explain why this enhancement would be useful** to many STAR users.
* **List some other tools where this enhancement exists.**

### Pull Requests

Please read the guides on creating good pull requests (PR) and follow these guidelines:
* **Use a clear and descriptive title** for the pull request.
* **State the purpose of the PR**: is it a bug-fix, new feature implementation, documentation improvement, or cosmetic change? Why is it important?
* **Explain the expected changes in STAR behavior**. Make sure that the default STAR behavior does not change.
* **Provide detailed code documentation and commit messages**.

Adopted from https://github.com/atom/atom/blob/master/CONTRIBUTING.md
