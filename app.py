from flask import Flask, render_template, request, jsonify
import re

app = Flask(__name__)

MAX_LEN = 30


def needleman_wunsch(seq_a, seq_b, match_score, mismatch_penalty, gap_penalty):
    n = len(seq_a)
    m = len(seq_b)

    # Initialize (n+1) x (m+1) DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    # Fill border gaps
    for i in range(n + 1):
        dp[i][0] = i * gap_penalty
    for j in range(m + 1):
        dp[0][j] = j * gap_penalty

    # Fill DP matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sc = match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_penalty
            dp[i][j] = max(
                dp[i - 1][j - 1] + sc,        # diagonal: match / mismatch
                dp[i - 1][j] + gap_penalty,    # up:       gap in seq_b
                dp[i][j - 1] + gap_penalty,    # left:     gap in seq_a
            )

    # Traceback from bottom-right to top-left
    path = []
    aligned_a, aligned_b = [], []
    i, j = n, m

    while i > 0 or j > 0:
        path.append([i, j])
        if i > 0 and j > 0:
            sc = match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_penalty
            if dp[i][j] == dp[i - 1][j - 1] + sc:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append(seq_b[j - 1])
                i -= 1
                j -= 1
            elif dp[i][j] == dp[i - 1][j] + gap_penalty:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append('-')
                i -= 1
            else:
                aligned_a.append('-')
                aligned_b.append(seq_b[j - 1])
                j -= 1
        elif i > 0:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append('-')
            i -= 1
        else:
            aligned_a.append('-')
            aligned_b.append(seq_b[j - 1])
            j -= 1

    path.append([0, 0])

    return (
        dp,
        path,
        ''.join(reversed(aligned_a)),
        ''.join(reversed(aligned_b)),
    )


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/align', methods=['POST'])
def align():
    data = request.get_json()

    seq_a = data.get('seq_a', '').strip().upper()
    seq_b = data.get('seq_b', '').strip().upper()

    try:
        match_score      = int(data.get('match', 1))
        mismatch_penalty = int(data.get('mismatch', -1))
        gap_penalty      = int(data.get('gap', -2))
    except (ValueError, TypeError):
        return jsonify({'error': 'Scores must be valid integers! Please check your inputs. 🔢'}), 400

    if not seq_a or not seq_b:
        return jsonify({'error': 'Please enter both sequences before aligning! 🧬'}), 400

    if len(seq_a) > MAX_LEN:
        return jsonify({
            'error': (f'Sequence A is {len(seq_a)} chars long — please keep both sequences '
                      f'under {MAX_LEN} characters for a clean visualization! ✂️')
        }), 400

    if len(seq_b) > MAX_LEN:
        return jsonify({
            'error': (f'Sequence B is {len(seq_b)} chars long — please keep both sequences '
                      f'under {MAX_LEN} characters for a clean visualization! ✂️')
        }), 400

    if not re.match(r'^[A-Z]+$', seq_a):
        return jsonify({'error': 'Sequence A should only contain letters A–Z (DNA, RNA & Amino Acids supported)! 🔬'}), 400

    if not re.match(r'^[A-Z]+$', seq_b):
        return jsonify({'error': 'Sequence B should only contain letters A–Z (DNA, RNA & Amino Acids supported)! 🔬'}), 400

    dp, path, aligned_a, aligned_b = needleman_wunsch(
        seq_a, seq_b, match_score, mismatch_penalty, gap_penalty
    )

    return jsonify({
        'matrix':    dp,
        'path':      path,
        'aligned_a': aligned_a,
        'aligned_b': aligned_b,
        'score':     dp[len(seq_a)][len(seq_b)],
        'seq_a':     seq_a,
        'seq_b':     seq_b,
    })


if __name__ == '__main__':
    app.run(debug=True)
