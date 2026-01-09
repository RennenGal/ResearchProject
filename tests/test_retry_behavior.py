"""
Property-based tests for retry behavior.

Tests the configurable retry mechanism with exponential backoff
to ensure it behaves correctly across all external database API calls.
"""

import pytest
import time
import asyncio
from unittest.mock import Mock, patch
from hypothesis import given, strategies as st, settings
from hypothesis.strategies import composite

from protein_data_collector.retry import RetryController, with_retry, with_retry_async
from protein_data_collector.config import RetryConfig


# Custom strategy for generating retry configurations
@composite
def retry_config_strategy(draw):
    """Generate valid RetryConfig instances for property testing."""
    max_retries = draw(st.integers(min_value=1, max_value=10))
    initial_delay = draw(st.floats(min_value=0.01, max_value=2.0))
    backoff_multiplier = draw(st.floats(min_value=1.1, max_value=5.0))
    max_delay = draw(st.floats(min_value=initial_delay, max_value=60.0))
    
    return RetryConfig(
        max_retries=max_retries,
        initial_delay=initial_delay,
        backoff_multiplier=backoff_multiplier,
        max_delay=max_delay
    )


class TestRetryBehavior:
    """Property-based tests for retry behavior."""
    
    @given(config=retry_config_strategy())
    @settings(max_examples=100, deadline=None)
    def test_configurable_retry_behavior_property(self, config):
        """
        Property 5: Configurable Retry Behavior
        
        For any external API call (InterPro, UniProt, MCP servers) that fails,
        the system should retry exactly K times (where K is configurable) using
        exponential backoff, then log the failure and continue processing.
        
        **Feature: protein-data-collector, Property 5: Configurable Retry Behavior**
        **Validates: Requirements 1.4, 3.8, 9.1, 9.2, 9.5**
        """
        controller = RetryController(config)
        
        # Track call attempts
        call_count = 0
        delays = []
        
        def failing_operation():
            nonlocal call_count
            call_count += 1
            raise ConnectionError(f"Simulated failure {call_count}")
        
        # Mock time.sleep to capture delays
        with patch('time.sleep') as mock_sleep:
            def capture_delay(delay):
                delays.append(delay)
            mock_sleep.side_effect = capture_delay
            
            # Execute operation that should fail after all retries
            with pytest.raises(ConnectionError):
                controller.execute_with_retry(
                    failing_operation,
                    database="TestDB",
                    operation_name="test_operation"
                )
        
        # Verify exactly K+1 attempts were made (initial + K retries)
        assert call_count == config.max_retries + 1
        
        # Verify exactly K delays were applied (no delay before first attempt)
        assert len(delays) == config.max_retries
        
        # Verify exponential backoff behavior
        for i, delay in enumerate(delays):
            expected_delay = config.initial_delay * (config.backoff_multiplier ** i)
            expected_delay = min(expected_delay, config.max_delay)
            assert abs(delay - expected_delay) < 0.001, f"Delay {i}: expected {expected_delay}, got {delay}"
        
        # Verify delays are non-decreasing (exponential backoff with max cap)
        for i in range(1, len(delays)):
            assert delays[i] >= delays[i-1] or delays[i] == config.max_delay
    
    @given(config=retry_config_strategy())
    @settings(max_examples=100, deadline=None)
    def test_successful_operation_no_retries(self, config):
        """
        Test that successful operations don't trigger retries.
        
        For any successful external API call, the system should return
        immediately without any retry attempts or delays.
        """
        controller = RetryController(config)
        
        call_count = 0
        expected_result = "success"
        
        def successful_operation():
            nonlocal call_count
            call_count += 1
            return expected_result
        
        # Mock time.sleep to ensure no delays occur
        with patch('time.sleep') as mock_sleep:
            result = controller.execute_with_retry(
                successful_operation,
                database="TestDB",
                operation_name="test_operation"
            )
        
        # Verify operation was called exactly once
        assert call_count == 1
        assert result == expected_result
        
        # Verify no delays were applied
        mock_sleep.assert_not_called()
    
    @given(
        config=retry_config_strategy(),
        success_attempt=st.integers(min_value=1, max_value=5)
    )
    @settings(max_examples=100, deadline=None)
    def test_eventual_success_within_retries(self, config, success_attempt):
        """
        Test that operations succeeding within retry limit work correctly.
        
        For any external API call that fails initially but succeeds within
        the retry limit, the system should return the successful result
        after the appropriate number of retry attempts.
        """
        # Ensure success_attempt is within the retry limit
        if success_attempt > config.max_retries + 1:
            success_attempt = config.max_retries + 1
        
        controller = RetryController(config)
        
        call_count = 0
        expected_result = "eventual_success"
        delays = []
        
        def eventually_successful_operation():
            nonlocal call_count
            call_count += 1
            if call_count < success_attempt:
                raise ConnectionError(f"Simulated failure {call_count}")
            return expected_result
        
        # Mock time.sleep to capture delays
        with patch('time.sleep') as mock_sleep:
            def capture_delay(delay):
                delays.append(delay)
            mock_sleep.side_effect = capture_delay
            
            result = controller.execute_with_retry(
                eventually_successful_operation,
                database="TestDB",
                operation_name="test_operation"
            )
        
        # Verify correct number of attempts
        assert call_count == success_attempt
        assert result == expected_result
        
        # Verify correct number of delays (success_attempt - 1)
        assert len(delays) == success_attempt - 1
        
        # Verify exponential backoff for the delays that occurred
        for i, delay in enumerate(delays):
            expected_delay = config.initial_delay * (config.backoff_multiplier ** i)
            expected_delay = min(expected_delay, config.max_delay)
            assert abs(delay - expected_delay) < 0.001
    
    @pytest.mark.asyncio
    @given(config=retry_config_strategy())
    @settings(max_examples=50, deadline=None)
    async def test_async_retry_behavior(self, config):
        """
        Test that async retry behavior matches synchronous behavior.
        
        For any async external API call that fails, the system should
        retry exactly K times using exponential backoff with async delays.
        """
        controller = RetryController(config)
        
        call_count = 0
        delays = []
        
        async def failing_async_operation():
            nonlocal call_count
            call_count += 1
            raise ConnectionError(f"Async failure {call_count}")
        
        # Mock asyncio.sleep to capture delays
        with patch('asyncio.sleep') as mock_sleep:
            async def capture_delay(delay):
                delays.append(delay)
            mock_sleep.side_effect = capture_delay
            
            with pytest.raises(ConnectionError):
                await controller.execute_with_retry_async(
                    failing_async_operation,
                    database="TestDB",
                    operation_name="async_test_operation"
                )
        
        # Verify exactly K+1 attempts were made
        assert call_count == config.max_retries + 1
        
        # Verify exactly K delays were applied
        assert len(delays) == config.max_retries
        
        # Verify exponential backoff behavior
        for i, delay in enumerate(delays):
            expected_delay = config.initial_delay * (config.backoff_multiplier ** i)
            expected_delay = min(expected_delay, config.max_delay)
            assert abs(delay - expected_delay) < 0.001
    
    @given(config=retry_config_strategy())
    @settings(max_examples=50, deadline=None)
    def test_decorator_retry_behavior(self, config):
        """
        Test that retry decorators behave consistently with direct controller usage.
        
        For any decorated function that fails, the retry decorator should
        apply the same retry logic as the controller.
        """
        call_count = 0
        delays = []
        
        @with_retry(database="TestDB", operation_name="decorated_operation", config=config)
        def decorated_failing_operation():
            nonlocal call_count
            call_count += 1
            raise ConnectionError(f"Decorated failure {call_count}")
        
        # Mock time.sleep to capture delays
        with patch('time.sleep') as mock_sleep:
            def capture_delay(delay):
                delays.append(delay)
            mock_sleep.side_effect = capture_delay
            
            with pytest.raises(ConnectionError):
                decorated_failing_operation()
        
        # Verify exactly K+1 attempts were made
        assert call_count == config.max_retries + 1
        
        # Verify exactly K delays were applied
        assert len(delays) == config.max_retries
        
        # Verify exponential backoff behavior
        for i, delay in enumerate(delays):
            expected_delay = config.initial_delay * (config.backoff_multiplier ** i)
            expected_delay = min(expected_delay, config.max_delay)
            assert abs(delay - expected_delay) < 0.001
    
    def test_delay_calculation_bounds(self):
        """
        Test that delay calculations respect max_delay bounds.
        
        For any retry configuration, calculated delays should never
        exceed the configured max_delay value.
        """
        config = RetryConfig(
            max_retries=10,
            initial_delay=1.0,
            backoff_multiplier=3.0,
            max_delay=5.0
        )
        
        controller = RetryController(config)
        
        # Test delay calculation for many attempts
        for attempt in range(15):
            delay = controller.calculate_delay(attempt)
            assert delay <= config.max_delay
            assert delay >= 0
        
        # Verify that without max_delay, delays would exceed the limit
        unbounded_delay = config.initial_delay * (config.backoff_multiplier ** 10)
        assert unbounded_delay > config.max_delay
        
        # But actual delay is capped
        actual_delay = controller.calculate_delay(10)
        assert actual_delay == config.max_delay