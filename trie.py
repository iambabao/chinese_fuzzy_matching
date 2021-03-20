# -*- coding: utf-8 -*-

"""
@Author             : sorahjy
@Date               : 2020/6/18 13:57
@Desc               :
@Last modified by   : Bao
@Last modified date : 2020/6/18 13:57
"""


def tokenize(seq):
    """
    change the tokenizer according to your language and application

    :param seq:
    :return:
    """

    return seq.split()


class Trie:
    def __init__(self):
        self.seqs = []
        self.seq2id = {}
        self.seq2tokens = {}

        # initialize with an empty root
        self.next_index = [{}]
        self.nodes = [[]]

    def add_seqs(self, seqs, skip=0):
        """

        :param seqs:
        :param skip:
        :return:
        """

        for seq in seqs:
            seq_id = len(self.seqs)
            tokens = tokenize(seq)

            self.seqs.append(seq)
            self.seq2id[seq] = seq_id
            self.seq2tokens[seq] = tokens
            for offset in range(min(skip + 1, len(tokens))):
                index = 0  # current depth
                for token in tokens[offset:]:
                    if token not in self.next_index[index]:  # add new node
                        self.next_index[index][token] = len(self.nodes)
                        self.next_index.append({})
                        self.nodes.append([])
                    index = self.next_index[index][token]
                self.nodes[index].append([seq_id, offset])

    def _approx(self, tokens, rep_pun=1.0, del_pun=1.0, add_pun=1.0, order_pun=0.00, max_pun=2.0, min_score=0.8):
        """

        :param tokens: input tokenized tokens
        :param rep_pun: punishment for replacement
        :param del_pun: punishment for deletion
        :param add_pun: punishment for addition
        :param order_pun: punishment for previous error
        :param max_pun: maximum punishment threshold
        :param min_score: minimum score threshold
        :return:
        """

        def _push(index, match, punishment):
            if visited.get((index, match), max_pun + 1e-6) > punishment:
                queue.append((index, match, punishment))
                visited[(index, match)] = punishment

        max_pun = min(max_pun, len(tokens) - 1)  # at least one token should be correct
        queue_head = 0
        queue = [(0, 0, 0)] # (index of tree nodes，current match of tokens，punishment of current match)
        matched_seqs = {}
        visited = {}

        first_token = True  # whether the first token should matched exactly
        while queue_head < len(queue):
            cur_index, cur_match, cur_pun = queue[queue_head]
            queue_head += 1
            if cur_match > len(tokens) or cur_pun > max_pun:
                continue
            if cur_match == len(tokens):
                for seq_id, offset in self.nodes[cur_index]:
                    seq = self.seqs[seq_id]
                    cur_pun += offset * del_pun
                    score = 1 - cur_pun / max(len(tokens), len(self.seq2tokens[seq]))
                    if score > min_score and score > matched_seqs.get(seq, 0):
                        matched_seqs[seq] = score
            cur_token = tokens[cur_match] if cur_match < len(tokens) else None  # move to next token
            next_index = self.next_index[cur_index].get(cur_token, -1)  # move to next node
            if next_index >= 0:
                # match token on tree
                _push(next_index, cur_match + 1, cur_pun)
            if not first_token:
                cur_order_pun = order_pun * max(len(tokens) - cur_match, 0)
                for token, next_index in self.next_index[cur_index].items():
                    # delete token on tree
                    _push(next_index, cur_match, cur_pun + del_pun + cur_order_pun)
                    # replace token on tree
                    if token != cur_token:
                        _push(next_index, cur_match + 1, cur_pun + rep_pun + cur_order_pun)
                # add token to tokens
                _push(cur_index, cur_match + 1, cur_pun + add_pun + cur_order_pun)
            first_token = False

        matched_seqs = sorted(matched_seqs.items(), key=lambda x: x[-1], reverse=True)
        matched_seqs = [(seq, score) for seq, score in matched_seqs]

        return matched_seqs

    def _exact(self, tokens):
        """

        :param tokens:
        :return:
        """

        def _push(index, match):
            if (index, match) not in visited:
                queue.append((index, match))
                visited.append((index, match))

        queue_head = 0
        queue = [(0, 0)] # (index of tree nodes，current match of tokens)
        matched_seqs = []
        visited = []

        while queue_head < len(queue):
            cur_index, cur_match = queue[queue_head]
            queue_head += 1
            if cur_match > len(tokens):
                continue
            if cur_match == len(tokens):
                for seq_id, offset in self.nodes[cur_index]:
                    if offset != 0:
                        continue
                    seq = self.seqs[seq_id]
                    matched_seqs.append(seq)
            for next_token, next_index in self.next_index[cur_index].items():
                if next_token == tokens[cur_match]:
                    _push(next_index, cur_match + 1)

        return matched_seqs

    def fuzzy_match(self, seq, **kwargs):
        """
        match input sequence with Trie Tree

        :param seq: untokenized input sequence
        :param kwargs:
        :return:
        """

        tokens = tokenize(seq)
        matched_seqs = self._approx(tokens, **kwargs)

        return matched_seqs

    def fuzzy_search(self, context, skip_overlap=True, **kwargs):
        """
        match each subsequence in context with Trie Tree

        :param context: untokenized input context
        :param skip_overlap: keep the longest match and skip overlapped subsequence if skip_overlap=True
        :param kwargs:
        :return:
        """

        tokens = tokenize(context)

        overlaps = 0
        matched_results = {}
        for i in range(len(tokens)):
            for j in range(len(tokens), i, -1):
                if skip_overlap and j <= overlaps:
                    break
                matched_seqs = self._approx(tokens[i:j], **kwargs)
                if len(matched_seqs) != 0:
                    matched_results[' '.join(tokens[i:j])] = matched_seqs
                    overlaps = j

        return matched_results


if __name__ == '__main__':
    trie = Trie()

    # add some commands
    trie.add_seqs([
        'import numpy',
        'import numpy as np',
        'from collections import Counter'
    ])

    # examples to check
    exams = [
        'import numpy',
        'import numpy as np',
        'import numpy as xx',
        'from collections import defaultdict',
    ]

    for e in exams:
        trie_match = trie.fuzzy_match(e, min_score=0.6)
        print('{}: {}'.format(e, trie_match))
    trie_match = trie.fuzzy_search('i first import numpy as np , then do my job.', min_score=0.6)
    print(trie_match)
    trie_match = trie.fuzzy_search('i first import numpy as np , then do my job.', skip_overlap=False, min_score=0.6)
    print(trie_match)
